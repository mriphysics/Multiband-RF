% Script to design singleband pulses.

clear all
close all
% Nt = 1024; %<-- Works well
Nt = 256;
tb = 8;

d1 = 0.01;
d2 = 0.01;

% % Linear-phase excitation - CVX not recommended.
% flip = pi/2;
% mode = 'cvx';
% type = 'exc';
% phase = 'linear';
% d1 = d1/10; %<-- Works.
 
% Linear-phase refocusing
flip = pi;
mode = 'cvx';
type = 'ref';
phase = 'linear';

% % Linear-phase matched-excitation refocusing pulse
% flip = pi;
% mode = 'cvx';
% type = 'me';
% phase = 'linear';

% % Minimum/maximum phase excitaion
% flip = pi/2;
% mode = 'cvx';
% type = 'exc';
% phase = 'minimum';

% % Minimum/maximum phase refocusing.
% flip = pi;
% mode = 'cvx';
% type = 'ref';
% phase = 'maximum';

% % Minimum/maximum phase matched-excitation refocusing.
% flip = pi;
% mode = 'cvx';
% type = 'me';
% phase = 'maximum';

% % quadratic phase excitation
% flip = pi/2;
% mode = 'cvx';
% type = 'exc';
% phase = 'quadratic';
% d1=d1/10;

% % quadratic phase refocusing
% flip = pi;
% mode = 'cvx';
% type = 'ref';
% phase = 'quadratic';

% % quadratic phase matched-excitation refocusing
% flip = pi;
% mode = 'cvx';
% type = 'me';
% phase = 'quadratic';

% % %  % Least-squres designs

% % Linear-phase excitation - Least-squares design
% flip = pi/2;
% mode = 'ls';
% type = 'exc';
% phase = 'linear';

% % Linear-phase refocusing 
% flip = pi;
% mode = 'ls';
% type = 'ref';
% phase = 'linear';

% % Linear-phase  matched-excitation 
% flip = pi;
% mode = 'ls';
% type = 'me';
% phase = 'linear';

% % minimum-phase excitation
% flip = pi/2;
% mode = 'ls';
% type = 'exc';
% phase = 'minimum';

% % minimum-phase refocusing
% flip = pi;
% mode = 'ls';
% type = 'ref';
% phase = 'minimum';

% % minimum-phase matched excitation
% flip = pi;
% mode = 'ls';
% type = 'me';
% phase = 'minimum';

% % quadratic-phase matched-excitation refocusing
% flip = pi;
% mode = 'ls';
% type = 'exc';
% phase = 'quadratic';

%  % quadratic-phase excitation
% flip = pi/2;
% mode = 'ls';
% type = 'exc';
% phase = 'quadratic';
 
% % quadratic-phase refocusing
% flip = pi;
% mode = 'ls';
% type = 'ref';
% phase = 'quadratic';
% % These parameters seem to work well!:
% d2 = d2/10; 
% Nt = 128;

% % quadratic-phase matched-excitation
% flip = pi;
% mode = 'ls';
% type = 'me';
% phase = 'quadratic';
% % These parameters seem to work well!:
% d2 = d2/10; 
% Nt = 120;

quiet = 1;

if strcmp(type,'me')
    [rf,tb,rf_me_exc] = singleband_rf(Nt,tb,flip,mode,type,phase,d1,d2,quiet);
else
    [rf,tb]           = singleband_rf(Nt,tb,flip,mode,type,phase,d1,d2,quiet);
end

slthick = 2*1e-3;
gam = 267522.1199722082;        % radians per sec per mT
b1max = 0.013;

dt = max(abs(rf))/(gam*b1max);
T = (length(rf))*dt;
t = (0:length(rf)-1)*dt;

BW = tb/T;
Gsel = 2*pi*BW/(gam*slthick)*ones(length(rf),1);

Nz = 10000;
xx = linspace(-0.02,0.02,Nz)';
pos = [0*xx 0*xx xx];
G =[0*Gsel 0*Gsel Gsel];

% if rescale_rf
%     rf=flip/sum(rf)*rf;
% end

switch type
    case 'exc'
        mxy_display = @(a,b)2*conj(a).*b;
        [~,~,~,~,as,bs] = blochsim_CK(rf(:)/gam/dt,G,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt);
        
        mxy = mxy_display(as(:,end),bs(:,end));
    case 'ref'
        mxy_display = @(a,b)b.^2;
        [~,~,~,~,as,bs] = blochsim_CK(rf(:)/gam/dt,G,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt);
        
        mxy = mxy_display(as(:,end),bs(:,end));
    case 'me'
        
%         % Method 1: Simulate as a single CK representation:
        mxy_display = @(a180,b180)(sqrt(2)*abs(b180).^4.*sqrt(1-abs(b180).^4/2));
        [~,~,~,~,as,bs] = blochsim_CK(rf(:)/gam/dt,G,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt);
        mxy = mxy_display(as(:,end),bs(:,end));
%         
% % %         Method 2: Simulate rf90 and rf180 separately.
%         dt90 = max(abs(rf_me_exc))/(gam*b1max);
%         T90 = (length(rf_me_exc))*dt90;        
% 
% %         BW90 = tb/T90; 
%         BW90 = 2*tb/T90; %<-- double BW corrects sb ripples. Why??
%         Gsel90 = 2*pi*BW90/(gam*slthick)*ones(length(rf_me_exc),1);
%         G90 =[0*Gsel90 0*Gsel90 Gsel90];
%         
%         [~,~,~,~,aexc,bexc] = blochsim_CK(rf_me_exc(:)/gam/dt90,G90,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt90);
%         [~,~,~,~,aref,bref] = blochsim_CK(rf(:)/gam/dt,G,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt);
%         
%         mxy = 2*conj(aexc(:,end)).*bexc(:,end).*bref(:,end).^2;
    otherwise        
        warning('pulse type not recognized - default to flip-angle representation');
        mxy_fa =@(a,b)acos(a.*conj(a)-b.*conj(b));
        mxy_display = mxy_fa;
end

% Allow option for flip_angle representation
if 0 
    mxy_fa =@(a,b)acos(a.*conj(a)-b.*conj(b));
    mxy_display = mxy_fa;
end

%%
fh = figure;
set(fh,'position',[243 415 1608 507]);
nr =1;
nc = 3;
% mxy = -bs(:,end).^2;

subplot(nr,nc,1)
hold on
grid on    
plot(t,real(rf)/dt/gam);
plot(t,imag(rf)/dt/gam);
ylabel('mT')
legend('Real B_{1,x}','Imag B_{1,y}');

subplot(nr,nc,2)
hold on;
plot(xx*1e3,abs(mxy))
hold on;
% Plot Passband constraints

% plot(xx*1e3,repmat([0.98:0.01:1],[length(xx) 1]),'--r','linewidth',0.5);
% Plot stopband constraints
% plot(xx*1e3,repmat(0:0.01:0.01,[length(xx) 1]),'--r','linewidth',0.5);

plot(xx*1e3,repmat([sin(flip/2)-2*d1, sin(flip/2)-d1, sin(flip/2)],[length(xx) 1]),'--r','linewidth',0.5);
plot(-slthick*1e3/2*ones(100,1),linspace(0, sin(flip/2),100),'--r','linewidth',0.5);
plot(+slthick*1e3/2*ones(100,1),linspace(0, sin(flip/2),100),'--r','linewidth',0.5);
grid on
xlabel('z/mm')
xlim([min(xx) max(xx)]*1e3)
xlim([-3*slthick*1e3 3*slthick*1e3]);

