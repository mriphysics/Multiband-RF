% Code to design either conventional or Amplitude-modulated (AM) only
% Mulitband RF pulses using 
% 1) Phase-optimzation (Wong, ISMRM proceedings 2012)
% 2) Time-shifting (Auerbach et al MRM 2013)
% 3) Root-flipping (Sharma et al MRM 2015)

% Requires CVX (http://cvxr.com/cvx/) and Pauly's RF tools
% (http://rsl.stanford.edu/research/software.html).

% Please use under MIT license (Copyright (c) 2015 mriphysics)
% Samy Abo Seada, Anthony Price, Jo Hajnal and Shaihan Malik. January 2017
% King's College London

clear all
close all;

gamma_mT = 267522.1; % radians/mT/s
b1max = 20*1e-3; % <--- peak B1 
mb = 6; % <-- Number of slices
tb = 4; % <-- Time-bandwidth product. Must be even number between 2 and 10.
bs = 10; % <-- Slice separation, in units of slice thicknesses
slthick = 2*1e-3;
AM_only = 1;    %<--- set to 0 for unconstrained design. 1 for AM design.

% Load in single-band refocusing pulses for matched-excitation
load('SB_SLR_cvxdesign_flip180_matchedexcitation.mat');

if or(mod(tb,2),or(tb<2,tb>10));
    fprintf('Time-bandwidth product must be even number between 2 and 10\n')
    fprintf('Alternatively, specify custom single-band waveform\n');
    keyboard;
end

% Obtain single-band waveform from the pulse structure. 
rfsb = pulse(tb/2).rf;
Nt = 2048; %<-- set number of time-points to interpolate to.
rfsb = length(rfsb)/Nt*interp1(linspace(0,1,length(rfsb)),rfsb,linspace(0,1,Nt));

% Spatial settings for Bloch-simulation.
FOV = 0.3; %<-- set Field-of-View in the slice-direction. in meters.
Nz = 4000; %<-- Number of spatial points.
z = linspace(-FOV/2,FOV/2,Nz)';
pos = [0*z 0*z z];

% Create autonomous functions with CK representation for refocusing and
% matched-excitation equations.

% For refocusings derivation see eq.9 in Pauly et al.  IEEE 1991
mxy_ref=@(b180)-b180.^2; 

% For matched-exciation expand eq.5 in Sharma et al. MRM 2015
mxy_me =@(b180)(sqrt(2)*abs(b180).^4.*sqrt(1-abs(b180).^4/2));

% By default, display matched-excitation profiles based on refocusing
% pulses.
mxy_display = mxy_me;

% Design Phase-optimized waveform.
rfmb_po = Phaseopt_fn(rfsb,mb,tb,bs,AM_only); %<-- in radians

% Find dwell-time such to scale to b1max.
dt_po = max(abs(rfmb_po))/gamma_mT/b1max; 

BW = tb/(length(rfmb_po)*dt_po);
Gsel = 2*pi*BW/(gamma_mT*slthick);

Gz = Gsel*ones(length(rfmb_po),1);
G = [0*Gz 0*Gz Gz];
% Cayley-Klein Bloch-simulation.
[~,~,~,~,~,b180_po] = blochsim_CK(rfmb_po/gamma_mT/dt_po,G,pos,ones(Nz,1),zeros(Nz,1),'dt',dt_po);
Mxy_po = mxy_display(b180_po(:,end));

% Design Time-shifted waveform
% Specify fractional shift. Note that for the paper results, we stepped
% through a range of fractional shifts (0:1/49:1) and reported the result
% for the time-shift for which the MB pulse had the minimum product between
% peak-amplitude and duration.
frac_shift = 0.5;
rfmb_ts = Timeshift_fn(rfsb,mb,tb,bs,frac_shift,AM_only); %<-- in radians.

% Find dwell-time such to scale to b1max.
dt_ts = max(abs(rfmb_ts))/gamma_mT/b1max;

% For time-shifting, the pulse duration is either the length of the MB
% pulse times the dwell-time, or single-band duration times
% (1+fractional-shift).
T_ts = length(rfmb_ts)*dt_ts;
T_ts = (1 + frac_shift)*length(rfsb)*dt_ts;

% Note that the BW for a time-shifted pulse is equivalent to its
% single-band waveform - see fig.1 in Auerbach et al. MRM 2013
BW = tb/(length(rfsb)*dt_ts);
Gsel = 2*pi*BW/(gamma_mT*slthick);

Gz = Gsel*ones(length(rfmb_ts),1);
G = [0*Gz 0*Gz Gz];
% Cayley-Klein Bloch-simulation.
[~,~,~,~,~,b180_ts] = blochsim_CK(rfmb_ts/gamma_mT/dt_ts,G,pos,ones(Nz,1),zeros(Nz,1),'dt',dt_ts);
Mxy_ts = mxy_display(b180_ts(:,end));

% Design Root-flipped waveform
N_rf = 512; %<-- Root-flipped design for many time-points can be time-consuming. Good idea to use less timepoints.
[rfmb_rf90,rfmb_rf180,tbsc]= rootflip_fn(N_rf,mb,tb,bs,AM_only);

dt_rf = max(abs(rfmb_rf180))/(gamma_mT*b1max);
BW = tbsc/(length(rfmb_rf180)*dt_rf);
Gsel = 2*pi*BW/(gamma_mT*slthick);        

Gz = Gsel*ones(length(rfmb_rf180),1);
G = [0*Gz 0*Gz Gz];

% Cayley-Klein Bloch-simulation.
[~,~,~,~,~,b180_rf] = blochsim_CK(rfmb_rf180/gamma_mT/dt_rf,G,pos,ones(Nz,1),zeros(Nz,1),'dt',dt_rf);
Mxy_rf = mxy_display(b180_rf(:,end));

%% Plot RF waveforms
% Find time-vectors [in ms] and convert multiband pulses from Radians to mT.
t_po = (0:length(rfmb_po)-1)*dt_po*1e3;
rfmb_po_mT = abs(rfmb_po)/gamma_mT/dt_po;

t_ts = (0:length(rfmb_ts)-1)*dt_ts*1e3;
rfmb_ts_mT = abs(rfmb_ts)/gamma_mT/dt_ts;

t_rf = (0:length(rfmb_rf180)-1)*dt_rf*1e3;
rfmb_rf_mT = abs(rfmb_rf180)/gamma_mT/dt_rf;

figure(1);
subplot(311); plot(t_po,rfmb_po_mT );
title(sprintf('Phase-optimizing T=%.2fms',t_po(end)));
xlabel('Time [ms]');ylabel('B_1 [mT]')
subplot(312); plot(t_ts,rfmb_ts_mT );
title(sprintf('Time-shifting T=%.2fms',t_ts(end)));
xlabel('Time [ms]');ylabel('B_1 [mT]')
subplot(313); plot(t_rf,rfmb_rf_mT );
title(sprintf('Root-flipping T=%.2fms',t_rf(end)));
xlabel('Time [ms]');ylabel('B_1 [mT]')

% Plot refocusing profile |Mxy|
figure(2);
plot(z,abs(Mxy_po));
hold on;grid on;
plot(z,abs(Mxy_ts));
plot(z,abs(Mxy_rf));
legend('Phase-Optimzing','Time-Shifting','Root-Flipping');

