clear all;
close all;

% Basic configurations - Load or design a new base SB pulse?
Nt = 2048; %<-- set number of time-points
gamma_mT = 2*pi*4.257*1e4; %<--- same as in minTime gradient function
load_SB_pulse = 0;  %<-- set to 1 to use precalculated Single-band pulse
flip = 180*pi/180; %<-- set flip-angle.

% ---- Set what type of multiband pulse. Choose from ---- %%
%    'no' : Non-optimized [Larkman 1991]
%    'po' : Phase-optimized [Hennig 1992, Wong 2013, Malik 2014]
%    'ts' : Time-shifting [ Auerbach MRM 2013]
%    'rf' : Root-flipping [Sharma MRM 2015]
%          The following four methods are described in Abo Seada MRM 2018
%    'mbv': Multiband (po) pulse followed by VERSE
%    'vmb': VERSE followed by multiband (Better for gradients!).
%    'mbvg' : Multiband VERSE pulse with GIRF-correction
%    'vmbg' : VERSE-multiband pulse with GIRF-correction 
%    'pins' : PINS pulses (Norris 2011)
%    'multipins' : MultiPINS pulses (Eichner 2014)

mb_type = 'vmbg'; 
mb = 4; %<-- Number of slices (Multiband factor)
tb = 6; %<-- Time-bandwidth product.
bs = 10;%<-- band separation (in units of slice-thicknesses)
slthick = 2*1e-3; %<-- slice thickness [mm]
gradientslopes = 1; %<-- Set to 1 for sloping gradients at start and end.

maxb1 = 13*1e-3; % <--- peak B1 [mT]
maxg=40; %<-- maximum gradient amplitude for VERSE [mT/m]
maxgslew = 200*1e3; %Maximum gradient slew-rate for VERSE [mT/m/s]
AM_only = 0; %<-- set to 1 for AM pulses

% % %     Select one of the three GIRFs    % % %
girf = load('h1_GIRF_20140729');disp('Using measured GIRF');
% girf = load('h2_GIRF_20170901.mat');disp('Using reconstructed GIRF');
% girf = [42*1e-6 42*1e-6];disp('Using analytical GIRF');
% % % ------------------------------------ % % %

return_gdem = 0; %<-- Options for GIRF-correction. 
                 %    If 1 the demanded gradient is returned.
if load_SB_pulse == 1
    % Example 1: Load pre-designed SB pulses
    load('SB_SLR_cvxdesign_flip180_quad_Mar27.mat');
    rfsb = pulse(tb-1).rf;                    

% %     % Example 2: Load in a vendor pulse
%     pulse = load('sg_XXX_XXX_X');
%     tb = 1;
%     rfsb = pulse.pulse;    
%     % Scale to FA
%     rfsb = FA*pi/180/sum(rfsb)*rfsb;
else
    % Example 2: Design a new pulse using the dz_singleband function
    Nt_dz = 300; %<-- Number of time-points used for single-band design
    d1 = 0.01;
    d2 = 0.01;
    
    mode = 'cvx'; %<-- set to cvx or ls
    pulse_type = 'me'; % Set to exc, ref or me (excitation, refocusing or matched-excitation refocusing respectively).
    phase = 'quadratic'; % Set to linear, minimum, maximum or quadratic.
    quiet=0;
    [rfsb,tb] = singleband_rf(Nt_dz,tb,flip,mode,pulse_type,phase,d1,d2,quiet);    
end

%% Use multiband_rf function
% Interploate to Nt
rfsb = length(rfsb)/Nt*interp1(linspace(0,1,length(rfsb)),rfsb,linspace(0,1,Nt))';                    

% Call main function
[rfmb,Gs,dtmb] = multiband_rf(mb_type,rfsb,mb,tb,bs,slthick,...
    maxb1,maxg,maxgslew,AM_only,girf,return_gdem,gradientslopes);

%% Run Bloch simulations.
plot_phase = 0; %<-- set to 1 to plot phase 
fprintf('\nRunning Bloch simulations...\n')
spos = (1:mb)-(mb+1)/2; 

Nz = 4096;
FOV = 6*(floor(mb/2)+1)*bs*slthick;
z = linspace(-FOV/2,FOV/2,Nz)';
pos = [z(:)*0 z(:)*0 z(:)];
[idz,~] = idmxy(z,slthick,mb,bs);
idz = idz(:);                        
mxy_exc =@(a,b)2*conj(a).*b;
mxy_fa =@(a,b)180/pi*acos(a.*conj(a)-b.*conj(b));
mxy_ref = @(a,b)b.^2;
if flip == pi
%         mxy_display = mxy_me;
    mxy_display = mxy_ref;
else                
    mxy_display = mxy_exc;
end

G3 = @(Gz) [0*Gz(:) 0*Gz(:) Gz(:)];

% Simulate Constant gradient MB pulse
[~,~,~,~,a,b] = blochsim_CK(rfmb,G3(Gs),pos,ones([Nz 1]),zeros([Nz 1]),'dt',dtmb);
mxy0 = mxy_display(a(:,end),b(:,end));

mxyzoom = [-2/3*mb*bs*slthick*1e2 2/3*mb*bs*slthick*1e2  -0.01 1.01];

t = (0:length(rfmb)-1)*dtmb*1e3;
fh = figure;
set(fh,'pos',[365 512 1287 459]);
subplot(2,2,1);plot(t,abs(rfmb)*1e3);
xlabel('Time [ms]');ylabel('B_1 [\mu T]');
title('RF pulse');

subplot(2,2,2);plot(t,Gs);
xlabel('Time [ms]');ylabel('Gs [mT/m]');
title('Gradient pulse');

subplot(2,2,[3 4]);plot(z*1e2,abs(mxy0));
xlabel('Position [cm]');ylabel('M_{xy} [a.u.]');
title('Slice profile');
axis(mxyzoom);

if plot_phase 
    demean = @(x) x(:)-mean(x(:))+0.5;
    [idz,pbidc] = idmxy(z,slthick,mb,bs);  
    for kk = 1:mb
        hold on;
        plot(z(pbidc(kk,:))*1e2,demean(angle(mxy0(pbidc(kk,:)))),'-','color',[0 0.5 0])        
    end
end
