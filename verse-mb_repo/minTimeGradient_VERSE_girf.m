function [C,time,g,slew,k, phi, sta, stb,p_of_t,VERSE_out] = minTimeGradient_VERSE_girf(C,g0, gfin, gmax, smax,T, ds, show,VERSE_in)
% 
% This function is based on minTimeGradientRIV part of the minTimeGradient
% package published by Michael Lustig, available here: 
%                       http://www.eecs.berkeley.edu/~mlustig/Software.html
%
% This modified version implements time optimal VERSE design published by
% Lee et al, MRM 2009 (doi: 10.1002/mrm.21950)
% 
% Coded by Shaihan Malik  (shaihan.malik@kcl.ac.uk) 2013, 2015
%
% [C,time,g,slew,k, phi, sta, stb,p_of_t,VERSE_out] =
%       minTimeGradient_VERSE(C,g0, gfin, gmax, smax,T, ds, show,VERSE_in)
%
% Inputs for the original function:
%
%   C       -   The Curve in k-space given in any parametrization [1/cm]
%               Accepts complex value C for 2D trajectory, and a Nx3 real
%               Matrix for 3D trajectories.
%   g0      -   Initial gradient amplitude (leave empty for g0 = 0)
%   gfin    -   Gradient value at the end of the trajectory. If not possible, 
%               the result would be the largest possible ampltude.
%               (Leave empty if you don't care to get maximum gradient.)
%   gmax    -   Maximum gradient [G/cm] (3.9 default)
%   smax    -   Maximum slew [G/Cm/ms]  (14.5 default)
%   T       -   Sampling time interval [ms] (4e-3 default)
%   ds      -   step size for ODE integration, leave empty to use default value
%   show    -   Show plots while optimizing (Warning: This will make the
%               process considerably slower!)
%   
% Return values from original function:
%   C       - reparametrized curve, sampled at T[ms]
%   time    - total time to get to the end
%   g       - gradiet waveform [G/cm]
%   s       - slew rate [G/cm/ms]
%   k       - exact k-space corresponding to gradient g (This function reparametrizes
%             C, then takes a derivative. Numerical errors in the derivative can lead to 
%             deviation.  
%   phi     - Geometry constraints on the amplitude vs. arclength
%   sta     - Solution for the forward ODE
%   stb     - Solution for the backward ODE
%
%   Additional input struct = VERSE_in
%                             VERSE_in.b = M x Nc RF pulse (Gauss)
%                             (M=#samples in time, Nc=#coils)
%                             VERSE_in.bmax = max B1 (Gauss)
%                             VERSE_in.Gmod = |G| (Gauss)
%                             VERSE_in.os = integer oversample factor (default 10)
%   
%   (oversample factor added to stabilize behaviour for complex RF pulses)
%
%   Additional output struct = VERSE_out
%                              VERSE_out.t = new time course sampled in specified interval (T)
%                              VERSE_out.bt= new pulse waveform(s) sampled in t
%                              VERSE_out.t1 = new time course oversampled by 
%                                        VERSE_in.os i.e. T/VERSE_in.os 
%                              VERSE_out.bt1= new pulse waveform(s) sampled in t1
% 

if nargin<1
    error('You gotta give me a curve!');
end

if nargin<2
    g0 = [];
end

if nargin<3
    gfin = [];
end

if nargin<4
    gmax = 4;
end

if nargin<5
    smax = 15;
end

if nargin<6
    T = 4e-3;
end

if nargin < 7
    ds = [];
end

if nargin<8
    show=0;
end

dt = T;

if isreal(C)
    if size(C,2) == 3
        C = complex3d(C);
    else
        error(' Curve can be complex Nx1 or real Nx3 matrix');
    end
end



gamma = 4.257;

% disp('         minTimeGradient_VERSE: Const arc-length parametrization');

% represent the curve using spline with parametrization p
Lp = length(C);
p = [0:Lp-1].';
PP = spline(p,C);


% interpolate curve for gradient accuracy
dp = 1e-1;
CC = ppval(PP, [0:dp:Lp-1].');


% find length of curve
Cp = (CC([2:end,end]) - CC)/dp;
Cp(end) = Cp(end-1);
Cp(1) = (CC(2)-CC(1))/dp;
s_of_p = cumtrapz(abs(Cp))*dp;
L = s_of_p(end);

%%% Shaihan: Interpolate pulse waveforms to get W
% disp('         minTimeGradient_VERSE: Extract RF pulse info for VERSE constraints');
b = interp1(p,VERSE_in.b,[0:dp:Lp-1]);
% make sure b is MxNc
if size(b,1)==1
    b = b.';
end
Gmod =  interp1(p,VERSE_in.Gmod,[0:dp:Lp-1]);

W = b./repmat(Gmod(:),[1 size(b,2)]);

bmax = VERSE_in.bmax;

% decide ds and compute st for the first point
stt0 = (gamma*smax) ; % always assumes first point is max slew
st0 = stt0*dt/2; % start at half the gradient for accuracy close to g=0
s0 = st0*dt;
if isempty(ds)
%     ds = s0/1.5; % smaller step size for numerical accuracy
    ds = s0/10; % sas MB VERSE smaller step size for numerical accuracy
end

s = [0 : ds : L].';
s_half = [0:ds/2:L].';
sta = s*0;

if isempty(g0)
    g0 = 0;
end
sta(1) = min(max(g0*gamma+st0), gamma*gmax);
p_of_s_half = interp1(s_of_p, [0:dp:Lp-1].', s_half, 'spline');
p_of_s = p_of_s_half(1:2:end);

%%% SJM interpolate W to s half
W_of_s_half = interp1([0:dp:Lp-1].',W, p_of_s_half); %<- again avoid spline
% disp('         minTimeGradient_VERSE: Compute geometry dependent constraints');
% compute constraints (forbidden line curve)
[phi,k] = sdotMax(PP,p_of_s_half,s_half, gmax, smax,W_of_s_half,bmax);
k = k([1:end,end,end]); % extend for the Runge-Kutte method


% disp('         minTimeGradient_VERSE: Solve ODE forward');
% solve ODE forward
for n=2:length(s)
    dstds = RungeKutte(s(n),ds,sta(n-1),k([(n-1)*2-1:(n-1)*2+1]),smax,L);
    tmpst = sta(n-1) + dstds;
    
    if isreal(tmpst) 
        sta(n) = min(tmpst,phi(n*2-1));
    else
        sta(n) = phi(n*2-1);
    end

    if mod(n,1000)==1 & show
       s(n)/L*100
       figure(100), plot(s,sta,linspace(s(1),s(end),length(phi)),phi); , axis([0,s(end),0,gmax*gamma*sqrt(2)]);, drawnow
   end
end


stb = 0*s;

if isempty(gfin)
    stb(end) = sta(end);
else
    stb(end) = min(max(gfin*gamma,st0), gamma*gmax);
end

% solve ODE backwards
% disp('         minTimeGradient_VERSE: Solve ODE backwards');
for n=length(s)-1:-1:1
    dstds = RungeKutte(s(n),ds,stb(n+1),k([(n+1)*2-1:-1:(n)*2-1]),smax,L);
    
    tmpst = stb(n+1) + dstds;
    
    if isreal(tmpst) 
        stb(n) = min(tmpst,phi(n*2-1));
    else
        stb(n) = phi(n*2-1);
    end

    if mod(n,1000)==1 & show
        s(n)/L*100
       figure(100), plot(s,stb,linspace(s(1),s(end), length(phi)),phi,s,sta);, axis([0,s(end),0,gmax*gamma*sqrt(2)]);,drawnow
   end
end

% disp('         minTimeGradient_VERSE: Final Interpolations')
% take the minimum of the curves 
st_of_s = min([sta,stb],[],2);

% compute time
t_of_s = cumtrapz(1./st_of_s*ds);
t = 0:dt:t_of_s(end);

s_of_t = interp1(t_of_s, s, t,'spline');
p_of_t = interp1(s, p_of_s, s_of_t,'spline');

C = ppval(PP,p_of_t);

if isa(C,'complex3d')
    C = get(C,'data');
else
    C = C(:);
end

g = diff(C)/gamma/dt; g=[g; g(end,:) + g(end,:)-g(end-1,:)];%(C([2:end, end]) - C([1:end-1, end-1]))/gamma/dt;
k = cumtrapz(g)*dt*gamma;

slew = diff(g)/dt;% (g([2:end, end]) - g([1:end-1, end-1]))/dt;
time = t(end);

% disp('         minTimeGradient_VERSE: Done');

%%% Shaihan: now re-interpolate W and hence make b at higher time resolution
if ~isfield(VERSE_in,'os')
    VERSE_in.os=10;
end
dt1 = dt/VERSE_in.os;
t1 = 0:dt1:t_of_s(end);
s_of_t1 = interp1(t_of_s, s, t1,'spline');
p_of_t1 = interp1(s, p_of_s, s_of_t1,'spline');
C1 = ppval(PP,p_of_t1);
C1 = get(C1,'data');
g1 = diff(C1)/gamma/dt1; g1=[g1; g1(end,:) + g1(end,:)-g1(end-1,:)];
g1abs = sum(g1.^2,2).^0.5;
% now interpolate W to this time base
W_of_t1 = interp1(s,W_of_s_half(1:2:end,:),s_of_t1);
if size(W_of_t1,1)==1%happens if Nc=1
    W_of_t1 = W_of_t1.';
end
bt1 = W_of_t1 .* repmat(g1abs,[1 size(W_of_t1,2)]);

%%% downsample to original time
% bt = interp1(t1,bt1,t);
% gt = interp1(t1,g1,t);

% sas 06Nov2016 investigate NaN cases.
bt = interp1(t1,bt1,t,'spline');
gt = interp1(t1,g1,t,'spline');

if any(isnan(bt))
    keyboard;
    warning('Found isnan in VERSE')
    bt(find(bt(isnan(bt))))=0;
end
VERSE_out.bt1=bt1;
VERSE_out.bt=bt;
VERSE_out.t=t;
VERSE_out.t1=t1;
VERSE_out.g=gt; %<-- sas added this ouput.



function [res] = RungeKutte(s,ds,st,k,smax,L)
% Solve ODE using Runge-Kutte method.

gamma = 4.257;
k = abs(k);
k1 = ds * 1/st*sqrt(gamma^2*smax^2 - abs(k(1))^2*st^4);
k2 = ds * 1/(st+ds*k1/2) * sqrt(gamma^2*smax^2 - abs(k(2))^2*(st+ds*k1/2)^4) ;
k3 = ds * 1/(st+ds*k2/2) * sqrt(gamma^2*smax^2 - abs(k(2))^2*(st+ds*k2/2)^4) ;
k4 = ds * 1/(st+ds*k3) * sqrt(gamma^2*smax^2 - abs(k(3))^2*(st+ds*k3)^4) ;

res = k1/6 + k2/3 + k3/3 + k4/6;



function [sdot, k] = sdotMax(PP, p_of_s, s, gmax, smax,W,bmax,alpha_int)

% [sdot, k, ] = sdotMax(PP, p_of_s, s, gmax, smax,W,bmax)
%
% Given a k-space curve C (in [1/cm] units), maximum gradient amplitude
% (in G/cm) and maximum slew-rate (in G/(cm*ms)).
% This function calculates the upper bound for the time parametrization
% sdot (which is a non scaled max gradient constaint) as a function of s.
% 
%   PP       --  spline polynomial
%   p_of_s  --  parametrization vs arclength
%   s       --  arclength parametrization (0->1)
%   gmax    --  maximum gradient (G/cm)
%   smax    --  maximum slew rate (G/ cm*ms)
%
%   returns the maximum sdot (1st derivative of s) as a function of arclength s
%   Also, returns curvature as a function of s and length of curve (L)
%
%  (c) Michael Lustig 2005 
%
%   Modified by Shaihan Malik 2015:
%
%   W       --  ratio b1/G
%   bmax    --  max B1 amplitude (Gauss)


gamma = 4.257;

s = s(:);
dp_p = p_of_s([2:end,end]) - p_of_s; , dp_p(end) = dp_p(end-1);
dp_m = p_of_s - p_of_s([1,1:end-1]);, dp_m(1) = dp_m(2);
ds_p = s([2:end,end]) - s; , ds_p(end) = ds_p(end-1);
ds_m = s - s([1,1:end-1]);, ds_m(1) = ds_m(2);


Cs_p = (ppval(PP,p_of_s + dp_p) - ppval(PP, p_of_s))./ds_p;
Cs_m = (ppval(PP,p_of_s) - ppval(PP, p_of_s-dp_m))./ds_m;
Cs = Cs_p/2 + Cs_m/2;
Css = (Cs_p - Cs_m)./(ds_m/2+ds_p/2);
k = abs(Css);
% fix edge numerical problems
k(end) = k(end-1);
k(1) = k(2);


% calc I constraint curve (maximum gradient)
% sdot1 = gamma*gmax*ones(size(s));

%%% Shaihan: Modify above, implement Eq.8 in Lee et al, 2009
sdot1a = gamma*gmax*ones(size(s));
%%% sas I think that to balance out the B1/W ratio
% sdot1a = gamma*gmax*ones(size(s))./alpha_int;
sdot1b = gamma*bmax./abs(W);%<- potentially multiple channels
% sas 07/10/2016 - add alpha for girf implementation
% sdot1b = gamma*bmax./abs(W)./alpha_int;%<- potentially multiple channels
sdot1 = min([sdot1a, sdot1b],[],2);

% calc II constraint curve (curve curvature dependent)
sdot2 = sqrt(gamma*smax ./ (abs(k)+eps));

% calc total constraint
sdot = min([sdot1, sdot2],[],2);

