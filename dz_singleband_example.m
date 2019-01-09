% Script to design singleband pulses.

clear all
close all

Nt = 256; %<-- nr of time-points
tb = 6;   %<-- Time bandwidth product

d1 = 0.01; %<-- passband ripple [%]
d2 = 0.01; %<-- stopband ripple [%]

% Set Flip angle, design mode, pulse type and phase type.
flip = pi; % Flip-angle in radians
mode = 'cvx'; % Set to cvx or ls
type = 'ref'; % Set to exc, ref or me (excitation, refocusing or matched-excitation refocusing respectively).
phase = 'linear'; % Set to linear, minimum, maximum or quadratic.
plot_fa = 0; %<-- set to 1 to plot profile in flip-angle representation
quiet = 1;  %<-- set to 1 to reduce command window output
slthick = 2*1e-3; %<-- Slice-thickness in m

if strcmp(type,'me')
    [rf,tb,rf_me_exc] = singleband_rf(Nt,tb,flip,mode,type,phase,d1,d2,quiet);
else
    [rf,tb]           = singleband_rf(Nt,tb,flip,mode,type,phase,d1,d2,quiet);
end

gamma_mT = 2*pi*4.257*1e4;
b1max = 0.013;

dt = max(abs(rf))/(gamma_mT*b1max);
T = (length(rf))*dt;
t = (0:length(rf)-1)*dt*1e3;

BW = tb/T;
Gs = 2*pi*BW/(gamma_mT*slthick)*ones(length(rf),1);

Nz = 10000;
xx = linspace(-0.02,0.02,Nz)';
pos = [0*xx 0*xx xx];
G =[0*Gs 0*Gs Gs];

%% Run Bloch simulations
fprintf('\nRunning Bloch simulations...\n')
switch type
    case 'exc'
        mxy_display = @(a,b)2*conj(a).*b;
        [~,~,~,~,as,bs] = blochsim_CK(rf(:)/gamma_mT/dt,G,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt);
        
        mxy = mxy_display(as(:,end),bs(:,end));
    case 'ref'
        mxy_display = @(a,b)b.^2;
        [~,~,~,~,as,bs] = blochsim_CK(rf(:)/gamma_mT/dt,G,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt);
        
        mxy = mxy_display(as(:,end),bs(:,end));
    case 'me'
        
%         % Method 1: Simulate as a single CK representation:
        mxy_display = @(a180,b180)(sqrt(2)*abs(b180).^4.*sqrt(1-abs(b180).^4/2));
        [~,~,~,~,as,bs] = blochsim_CK(rf(:)/gamma_mT/dt,G,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt);
        mxy = mxy_display(as(:,end),bs(:,end));
         
% % %         Method 2: Simulate rf90 and rf180 separately.
%         dt90 = max(abs(rf_me_exc))/(gamma_mT*b1max);
%         T90 = (length(rf_me_exc))*dt90;        
% 
% %         BW90 = tb/T90; 
%         BW90 = 2*tb/T90; %<-- double BW corrects sb ripples. Why??
%         Gsel90 = 2*pi*BW90/(gamma_mT*slthick)*ones(length(rf_me_exc),1);
%         G90 =[0*Gsel90 0*Gsel90 Gsel90];
%         
%         [~,~,~,~,aexc,bexc] = blochsim_CK(rf_me_exc(:)/gamma_mT/dt90,G90,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt90);
%         [~,~,~,~,aref,bref] = blochsim_CK(rf(:)/gamma_mT/dt,G,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt);
%         
%         mxy = 2*conj(aexc(:,end)).*bexc(:,end).*bref(:,end).^2;
    otherwise        
        warning('pulse type not recognized - default to flip-angle representation');
        mxy_fa =@(a,b)acos(a.*conj(a)-b.*conj(b));
        mxy_display = mxy_fa;
end

% Allow option for flip_angle representation
if plot_fa
    mxy_fa =@(a,b)acos(a.*conj(a)-b.*conj(b));
    mxy_display = mxy_fa;
    mxy = mxy_display(as(:,end),bs(:,end));
end

%%
fh = figure;
nr =2;
nc =2;

subplot(nr,nc,1)
hold on
grid on    
plot(t,real(rf)/dt/gamma_mT*1e3);
plot(t,imag(rf)/dt/gamma_mT*1e3);
xlabel('Time [ms]');ylabel('B_1 [\mu T]');
title('RF pulse');
legend('Real','Imag');

subplot(nr,nc,2)
hold on
grid on    
plot(t,Gs);
xlabel('Time [ms]');ylabel('Gs [mT/m]');
title('Gradient pulse');

subplot(nr,nc,[3 4])
hold on;
plot(xx*1e3,abs(mxy))
% hold on;
% plot(xx*1e3,repmat([sin(flip/2)-2*d1, sin(flip/2)-d1, sin(flip/2)],[length(xx) 1]),'--r','linewidth',0.5);
% plot(-slthick*1e3/2*ones(100,1),linspace(0, sin(flip/2),100),'--r','linewidth',0.5);
% plot(+slthick*1e3/2*ones(100,1),linspace(0, sin(flip/2),100),'--r','linewidth',0.5);
xlabel('Position [cm]');ylabel('M_{xy} [a.u.]');
title('Slice profile');
xlim([min(xx) max(xx)]*1e3)
xlim([-3*slthick*1e3 3*slthick*1e3]);

%% Testbench: Design settings for number of cases

% % Linear-phase excitation - CVX not recommended.
% flip = pi/2;
% mode = 'cvx';
% type = 'exc';
% phase = 'linear';
% d1 = d1/10; %<-- Works.
%  
% % Linear-phase refocusing
% flip = pi;
% mode = 'cvx';
% type = 'ref';
% phase = 'linear';
% 
% % Linear-phase matched-excitation refocusing pulse
% flip = pi;
% mode = 'cvx';
% type = 'me';
% phase = 'linear';
% 
% % Minimum/maximum phase excitaion
% flip = pi/2;
% mode = 'cvx';
% type = 'exc';
% phase = 'minimum';
% 
% % Minimum/maximum phase refocusing.
% flip = pi;
% mode = 'cvx';
% type = 'ref';
% phase = 'maximum';
% 
% % Minimum/maximum phase matched-excitation refocusing.
% flip = pi;
% mode = 'cvx';
% type = 'me';
% phase = 'maximum';
% 
% % quadratic phase excitation
% flip = pi/2;
% mode = 'cvx';
% type = 'exc';
% phase = 'quadratic';
% d1=d1/10;
% 
% % quadratic phase refocusing
% flip = pi;
% mode = 'cvx';
% type = 'ref';
% phase = 'quadratic';
% 
% % quadratic phase matched-excitation refocusing
% flip = pi;
% mode = 'cvx';
% type = 'me';
% phase = 'quadratic';
% 
% % %  % Least-squres designs
% 
% % Linear-phase excitation - Least-squares design
% flip = pi/2;
% mode = 'ls';
% type = 'exc';
% phase = 'linear';
% 
% % Linear-phase refocusing 
% flip = pi;
% mode = 'ls';
% type = 'ref';
% phase = 'linear';
% 
% % Linear-phase  matched-excitation 
% flip = pi;
% mode = 'ls';
% type = 'me';
% phase = 'linear';
% 
% % minimum-phase excitation
% flip = pi/2;
% mode = 'ls';
% type = 'exc';
% phase = 'minimum';
% 
% % minimum-phase refocusing
% flip = pi;
% mode = 'ls';
% type = 'ref';
% phase = 'minimum';
% 
% % minimum-phase matched excitation
% flip = pi;
% mode = 'ls';
% type = 'me';
% phase = 'minimum';
% 
% % quadratic-phase matched-excitation refocusing
% flip = pi;
% mode = 'ls';
% type = 'exc';
% phase = 'quadratic';
% 
%  % quadratic-phase excitation
% flip = pi/2;
% mode = 'ls';
% type = 'exc';
% phase = 'quadratic';
%  
% % quadratic-phase refocusing
% flip = pi;
% mode = 'ls';
% type = 'ref';
% phase = 'quadratic';
% % These parameters seem to work well!:
% d2 = d2/10; 
% Nt = 128;
% 
% % quadratic-phase matched-excitation
% flip = pi;
% mode = 'ls';
% type = 'me';
% phase = 'quadratic';
% % These parameters seem to work well!:
% d2 = d2/10; 
% Nt = 120;
