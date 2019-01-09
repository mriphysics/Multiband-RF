clear all;
close all;

% 23/04/2018 - sas
% This script will design five different Multiband pulses (Constant
% gradient, vMB, MBv, PINS and MultiPINS), and simulate their expected
% gradient distortion based on a selected Gradient Impulse Response
% Function (GIRF).

% Use this top section to modify any design parameters and select the GIRF
% you would like to use.

% This code has been released in combination with our publication in MRM,
% which is released with a CC-BY license.
% This mean: feel free to modify and spread as long as acknowledgement towards
% the original authors is made.

% Samy Abo Seada
% 2018 King's College London
Nt = 2048;                 %<--- Nr of time-points
gamma_mT = 2*pi*4.257*1e4; %<--- Gyromagnetric ratio [rad/mT/s]
slthick = 2*1e-3;          %<--- Slice thickness [m]
b1max = 13*1e-3;           %<--- Peak B1 amplitude [mT]
maxg = 40;                 %<--- Maximum gradient amplitude [mT/m]
maxgslew = 200*1e3;        %<--- Maximum gradient amplitude [mT/m/s]
AM_only = 0;               %<--- Enable AM-only

% --- Pulse characteristics --- %
mb = 3;                    %<--- Multiband factor (nr of slices)
tb = 4;                    %<--- time bandwidth product (integer from 2 to 10)
bs = 14;                   %<--- Slice gap in units of slice-thickness [dimensionless].

% -- Select whether to plot refocusing profiles (1) or in pure Flip-angle (0) --%
plot_refocusing = 0; 

% --- Select GIRF --- %
girf_idc = 1;     %<--- set to 1 for a measured Philips GIRF
                        %      2 for a reconstructed Siemens GIRF
                        %      3 for analytical Lorentzian-shape GIRF
%% Specify GIRF
switch girf_idc
    case 1
    % h1_GIRF_20140729 is a measured Philips GIRF used for experimental results.
    girf = load('bin/h1_GIRF_20140729');
    structural_girf = 1;
    case 2
    % h2_GIRF_20170901 is a reconstructed GIRF from Testud, using an
    % Siemens gradient system with a higher temporal BW.
    girf = load('bin/h2_GIRF_20170901');
    structural_girf = 1;
    case 3

    % It's possible to specify an analytical GIRF model. Here, a
    % mono-exponential filter will roughly approximate h1 to have the frequency
    % response of a Lorentzian, with a fixed delay-term.

    tau = 41.717*1e-6; %<-- Time-constant for mono-exponential filter
    delay =  45.442*1e-6; %<-- constant delay.
    girf = [tau delay];
    structural_girf = 0;
end

%% Design Singleband pulse
load('SB_SLR_cvxdesign_flip180_quad_Mar27.mat')
rfsb = pulse(tb-1).rf;                    

rfsb = length(rfsb)/Nt*interp1(linspace(0,1,length(rfsb)),rfsb,linspace(0,1,Nt))';                    

dt_sb = max(abs(rfsb))/(gamma_mT*b1max);
rfsb_mT = rfsb ./(gamma_mT*dt_sb);

BW_sb = tb/(length(rfsb_mT)*dt_sb);
Gsel_sb = 2*pi*BW_sb/(gamma_mT*slthick);                                      
Gz_sb = Gsel_sb*ones(length(rfsb_mT),1);

%% Design Phase-optimized MB pulse
rfmb = Phaseopt_fn(rfsb,mb,tb,bs,AM_only);

dt=max(abs(rfmb))./(gamma_mT*b1max);
% Convert MB pulse from rad to mT:
rfmb = rfmb./(gamma_mT*dt);
BW = tb/(Nt*dt);
Gsel = 2*pi*BW/(gamma_mT*slthick);
Gz = Gsel*ones(Nt,1);
gmb = [0*Gz 0*Gz Gz];    
t = (0:length(rfmb)-1)*dt;


%% Design linear-phase MBv pulse

fprintf('Designing MBv pulse...\n')
verse_singleband = 0;
dt_os = 3;

[rfMBv,gMBv,gMBv_actual]= dz_MBverse(rfmb,Gz,dt,maxg,...
    maxgslew,b1max,verse_singleband,mb,bs*slthick,dt_os,AM_only,girf);
Nv = length(rfMBv);
dtv =dt/dt_os;
tv = 0:dtv:(Nv-1)*dtv;

%% Design linear-phase vMB pulse
fprintf('Designing vMB pulse...\n')

rf_init = rfsb_mT;

% Scale input single-band RF to mT:
verse_singleband = 1;
dt_os = 2;
epsilon = 1e-4;
max_iterations = 1;

[rfvMB,gvMB,gvMB_actual]= dz_MBverse(rf_init,Gz_sb,dt_sb,maxg,...
    maxgslew,b1max,verse_singleband,mb,bs*slthick,dt_os,AM_only,girf);
Nvmb = length(rfvMB);
dtvmb =dt_sb/dt_os;
tvmb = 0:dtvmb:(Nvmb-1)*dtvmb;

%% Design PINS pulse

fprintf('Designing PINS pulse...\n')
halfShift = ~mod(mb,2); % shift pattern by 1/2 slice gap to line up with target
mindurRF = 1;
de = 0.01;

[pins_rf,pins_g] = dz_pins_arbSB(gamma_mT*dt_sb*rfsb_mT,tb,bs*slthick,slthick,de,b1max,maxg,maxgslew,...
     dt_sb,mindurRF,halfShift);

t_pins = (0:length(pins_rf)-1)*dt_sb;
pins_G = [0*pins_g 0*pins_g pins_g];

if structural_girf
    pins_Gactual = gradient_distort_GIRF(pins_G,girf.ff,girf.Hw,dt_sb,500);
else %<-- Use mono-exponential GIRF
    pins_Gactual = gradient_distort_FT(pins_G,girf(1),girf(2)*ones(1,3),dt_sb,500);
end

%% Design Time-optimal MultiPINS pulse
% Input to Time_Optimal_Multipins is single-band pulse scaled in radians.
fprintf('Designing MultiPINS pulse...\n')
halfShift = ~mod(mb,2); % shift pattern by 1/2 slice gap to line up with target
mindurRF = 1; % switch to use min duration RF for all subpulses
Mixing_ratio_v = 0:0.005:1;

[ Mupins_rf,Mupins_g,Mixing_ratio,RF_energy_vs_mixingRatio] = Time_Optimal_Multipins(...
    gamma_mT*dt_sb*rfsb_mT,mb,tb,bs*slthick,slthick,de,b1max,maxg,maxgslew,dt_sb,mindurRF,halfShift,Mixing_ratio_v);

Mupins_G = [0*Mupins_g 0*Mupins_g Mupins_g];

if structural_girf
    Mupins_Gactual = gradient_distort_GIRF(Mupins_G,girf.ff,girf.Hw,dt_sb,500);
else %<-- Use mono-exponential GIRF
    Mupins_Gactual = gradient_distort_FT(Mupins_G,girf(1),girf(2)*ones(1,3),dt_sb,500);
end

t_Mupins = (0:length(Mupins_rf)-1)*dt_sb;       

%% Run Bloch simulations.
fprintf('\nRunning Bloch simulations...\n')
spos = (1:mb)-(mb+1)/2; 

Nz = 4096;
FOV = 3*(floor(mb/2)+1)*bs*slthick;
z = linspace(-FOV/2,FOV/2,Nz)';

pos = [z(:)*0 z(:)*0 z(:)];

mxy_fa =@(a,b)180/pi*acos(a.*conj(a)-b.*conj(b));
mxy_ref = @(a,b)b.^2;

if plot_refocusing
    mxy_display = mxy_ref;
else
    mxy_display = mxy_fa;
end

% Simulate Constant gradient MB pulse
fprintf('Running for constant gradient MB pulses...\n')
[~,~,~,~,a,b] = blochsim_CK(rfmb,gmb,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt);
mxy0 = mxy_display(a(:,end),b(:,end));

fprintf('Running for MBv pulses ...\n')
[~,~,~,~,a,b] = blochsim_CK(rfMBv,gMBv,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dtv);
mxy_mbv = mxy_display(a(:,end),b(:,end));
[~,~,~,~,a,b] = blochsim_CK(rfMBv,gMBv_actual,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dtv);
mxy_mbv_act = mxy_display(a(:,end),b(:,end));

fprintf('Running for vMB pulses...\n')
[~,~,~,~,a,b] = blochsim_CK(rfvMB,gvMB,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dtvmb);
mxy_vmb = mxy_display(a(:,end),b(:,end));
[~,~,~,~,a,b] = blochsim_CK(rfvMB,gvMB_actual,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dtvmb);
mxy_vmb_act = mxy_display(a(:,end),b(:,end));

fprintf('Running for PINS pulses...\n')
[~,~,~,~,a,b] = blochsim_CK(pins_rf,pins_G,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt_sb);
mxy_pins = mxy_display(a(:,end),b(:,end));
[~,~,~,~,a,b] = blochsim_CK(pins_rf,pins_Gactual,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt_sb);
mxy_pins_act = mxy_display(a(:,end),b(:,end));

fprintf('Running for MultiPINS pulses...\n')
[~,~,~,~,a,b] = blochsim_CK(Mupins_rf,Mupins_G,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt_sb);
mxy_multipins = mxy_display(a(:,end),b(:,end));
[~,~,~,~,a,b] = blochsim_CK(Mupins_rf,Mupins_Gactual,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt_sb);
mxy_multipins_act = mxy_display(a(:,end),b(:,end));

%% Plot RF, gradient waveforms and slice-profiles
fh = figure;
nr = 5;
nc = 5;

% Plot constant-gradient MB
subplot(nr,nc,0*nc+1);plot(t*1e3,abs(rfmb)*1e3);
ylabel('MB (const)');

subplot(nr,nc,0*nc+2);plot(t*1e3,abs(gmb(:,3)));
subplot(nr,nc,0*nc+3:nc);plot(z*1e2,abs(mxy0));

% Plot MB verse
subplot(nr,nc,1*nc+1);plot(tv*1e3,abs(rfMBv)*1e3);
ylabel('MBv');

subplot(nr,nc,1*nc+2);plot(tv*1e3,abs(gMBv(:,3)));
hold on;plot(tv*1e3,abs(gMBv_actual(:,3)));

subplot(nr,nc,1*nc+3:2*nc);plot(z*1e2,abs(mxy_mbv));
hold on;plot(z*1e2,abs(mxy_mbv_act));
legend('Target','Predicted gradient distortion');

% Plot verse MB
subplot(nr,nc,2*nc+1);plot(tvmb*1e3,abs(rfvMB)*1e3);
ylabel('vMB');

subplot(nr,nc,2*nc+2);plot(tvmb*1e3,abs(gvMB(:,3)));
hold on;plot(tvmb*1e3,abs(gvMB_actual(:,3)));

subplot(nr,nc,2*nc+3:3*nc);plot(z*1e2,abs(mxy_vmb));
hold on;plot(z*1e2,abs(mxy_vmb_act));

% Plot PINS
subplot(nr,nc,3*nc+1);plot(t_pins*1e3,abs(pins_rf)*1e3);
ylabel('PINS');

subplot(nr,nc,3*nc+2);plot(t_pins*1e3,abs(pins_G(:,3)));
hold on;plot(t_pins*1e3,abs(pins_Gactual(:,3)));

subplot(nr,nc,3*nc+3:4*nc);plot(z*1e2,abs(mxy_pins));
hold on;plot(z*1e2,abs(mxy_pins_act));

% Plot multiPINS
subplot(nr,nc,4*nc+1);plot(t_Mupins*1e3,abs(Mupins_rf)*1e3);
ylabel('MBv');

subplot(nr,nc,4*nc+2);plot(t_Mupins*1e3,abs(Mupins_G(:,3)));
hold on;plot(t_Mupins*1e3,abs(Mupins_Gactual(:,3)));

subplot(nr,nc,4*nc+3:5*nc);plot(z*1e2,abs(mxy_multipins));
hold on;plot(z*1e2,abs(mxy_multipins_act));
