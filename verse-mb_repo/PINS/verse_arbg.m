function [ RFout ] = verse_arbg( RFin, Gin, dtin, Garb,dtarb,smax)

% 19/10/16 Samy Abo Seada - Reshapes RF waveform to preserve the resultant slice profile
% for a combination of RF and gradient input, given an arbitrary gradient
% Garb. 
% This can even work when Garb is defined on a different time-base.

% This is done by simultaneously preserving the incremental change in
% k-space, and preserving W(s), which is the ratio arclength change of B1/G at every
% time-step.

% Initially developed in script MB_verse_mappingscript for work on
% Multiband VERSE GIRF.

% 23/04/2018 sas - version to release as part of verse-mb publication. 

if nargin <6
    smax = 2e5; %<-- gradient slew rate in mT/m/s
end
gamma_mT = 2*pi*4.257*1e4; %<--- same as in minTime gradient function

Nt = length(RFin);
t = dtin*(1:Nt)';

Ntarb = length(Garb);
tarb = dtarb*(1:Ntarb)';

kin = gamma_mT*cumtrapz(t,Gin);

p = 0:1:Nt-1; %<--- path on normalised sampled grid
PP = spline(p,kin);

% interpolate curve for gradient accuracy
dp = 1e-1;
ph = 0:dp:Nt-1; %<--- path on highly sampled grid. Use for interpolating RFin, Gin and W.
CC = ppval(PP, ph'); %<--- piece-wise polynomial fitting

% find length of curve
Cp = (CC([2:end,end]) - CC)/dp;
Cp(end) = Cp(end-1);
Cp(1) = (CC(2)-CC(1))/dp;
s_of_p = cumtrapz(abs(Cp))*dp;
L = s_of_p(end);

% Interpolate input waveforms onto the oversampled path
B1_of_p = interp1(p,RFin,ph);
G_of_p = interp1(p,Gin,ph);
W_of_p = B1_of_p./G_of_p;

% Convert W(p) to W(s)
ds = gamma_mT*smax*dtin/2*dtin/10; %
s = [0:ds:L].';
p_of_s = interp1(s_of_p,ph,s);
W_of_s = interp1(ph,W_of_p,p_of_s);

% The next bit of code reverses the code in the TO-VERSE algorithm, by
% going from G(t) to s_of_t.
% s_of_t is just the integral of the output gradient.

s_of_t = gamma_mT*cumtrapz(tarb,Garb(:));
W_of_t = interp1(s,W_of_s,s_of_t);

% 27/02/2017 - Spline gives numerical errors in my pulses. linear
% interpolation leads to NaN at the final 1-5 elements after the
% undersampling from s to s(t). Replace NaN with zero.
W_of_t(isnan(W_of_t)) = 0;
RFout = W_of_t.*repmat(Garb(:),[1 size(W_of_t,2)]);

end

