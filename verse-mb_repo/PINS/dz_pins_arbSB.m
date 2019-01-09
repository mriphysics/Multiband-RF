function [rf,g] = dz_pins_arbSB(rfsb,tb,slsep,slthick,de,maxb1,gmax,gslew,dt,minRFdur,halfShift)

% design a PINS refocusing pulse
% 10/1/12 WA Grissom

% Output variables are: rf, g

% 26/08/2016 sas - Adapt Grissoms fn for arbitrary soft-pulse specification and use SI
% units.
% Version 2 allows for arbitrary sb pulse input.
% Original source: https://bitbucket.org/wgrissom/lowpeakpowermbrf

% 23/04/2018 sas - version to release as part of verse-mb publication. 
if ~exist('halfShift','var')
  halfShift = false; % do not shift profile by half a slice gap
end

d1 = de/4; % ripple in ref passband profile, converted to beta
d2 = sqrt(de); % ripple in ref stopband profile, converted to beta

kzw = tb/slthick; % 1/m, width in kz-space we must go
Npulses = ceil(kzw/(1/slsep)); % number of subpulses

% call SLR to get pulse
% rfSoft = imag(dzrf(Npulses,tb,'se','ls',d1,d2));
% cplot(rfSoft);
% 28/08/16 Allow arbitrary sbpulse input. The following line downsamples
% the specified pulse to the number of PINS subpulses
rfSoft = interp1(linspace(0,1,length(rfsb)),rfsb,linspace(0,1,Npulses));

% sas 21/02/17 scale back to original flip angle.
rfSoft = rfSoft*sum(rfsb)/sum(rfSoft);

if halfShift
  rfSoft = rfSoft.*((-1).^[0:Npulses-1]);
end

gam = 2*pi*4.257*1e4; %<--- same as in minTime gradient function

% design the blip trapezoid. Units irrelevant for dotrap fn.
garea = 1/slsep/(gam/2/pi); % mT/m*s, k-area of each blip
gzblip = dotrap(garea,gmax,gslew,dt).'; 

if ~minRFdur

    % matched-duration RF subpulses
%     hpw = ceil(max(abs(rfSoft))./(2*pi*gambar*maxb1*dt));
    hpw = ceil(max(abs(rfSoft))./(gam*maxb1*dt));

    % interleave RF subpulses with gradient blips to form full pulses
    rf = kron(rfSoft(:),[ones(hpw,1);zeros(size(gzblip))]);

    g = repmat([zeros(hpw,1);gzblip(:)],[Npulses 1]);
    
    rf = rf./(sum(rf(:))*gam*dt)*sum(rfSoft); % convert to mT

else

    % matched-amplitude RF subpulses
    rf = [];g = [];
    for ii = 1:Npulses

        % Find non-integer (ni) number of subpulse-samples
        hpw_ni = abs(rfSoft(ii))/(gam*maxb1*dt);
        hpw = ceil(abs(rfSoft(ii))./(gam*maxb1*dt));
        
        beta =maxb1/hpw * (hpw-hpw_ni);
        rf = [rf; sign(rfSoft(ii))*gam*dt*(maxb1-beta)*ones(hpw,1);zeros(size(gzblip))];
        g = [g; zeros(hpw,1); gzblip(:)];
    end   
    
    rf = rf./(sum(rf(:))*gam*dt)*sum(rfSoft); % convert to mT
end

% remove the last blip, if 'se'
g = g(1:end-length(gzblip)); % flip sign (gradient reversal to suppress fat) and remove last blip
rf = rf(1:end-length(gzblip));
