function [rf,g] = dz_Multipins(rfsb,mb,tb,slsep,slthick,de,maxb1,gmax,gslew,dt,minRFdur,halfShift,Mixing_ratio)

% design a PINS refocusing pulse
% 10/1/12 WA Grissom

% Output variables are: rf, g

% 26/08/2016 sas Use Grissoms code as a basis to design MultiPINS pulses
% (Eichner2014). 
% Original source: https://bitbucket.org/wgrissom/lowpeakpowermbrf

% 23/04/2018 sas - version to release as part of verse-mb publication. 

if ~exist('halfShift','var')
  halfShift = false; % do not shift profile by half a slice gap
end

kzw = tb/slthick; % 1/m, width in kz-space we must go
Npulses = ceil(kzw/(1/slsep)); % number of subpulses

if isempty(rfsb)
    d1 = de/4; % ripple in ref passband profile, converted to beta
    d2 = sqrt(de); % ripple in ref stopband profile, converted to beta

    %  call SLR to get pulse
    rfSoft = imag(dzrf(Npulses,tb,'se','ls',d1,d2));
    % cplot(rfSoft);
else% 28/08/16 Allow arbitrary sbpulse input. The following line downsamples
    % the specified pulse to the number of PINS subpulses
    rfSoft = length(rfsb)/Npulses*interp1(linspace(0,1,length(rfsb)),rfsb,linspace(0,1,Npulses));
    
    % sas 21/02/17 scale back to original flip angle. Found that this is
    % slightly better for Refocusing pulses, but works evenly well without
    % rescaling.
    rfSoft = rfSoft*sum(rfsb)/sum(rfSoft);
end

%% sas 26/08/16 Create a Multiband pulse

gam = 2*pi*4.257*1e4; %<--- same as in minTime gradient function

BW = tb/(length(rfsb)*dt);
% Find Gsel
Gsel = 2*pi*BW./(gam*slthick);
% Find modulation function
spos = (1:mb)-(mb+1)/2; 
t = (1:length(rfsb))*dt;t=t(:);
phi_sel = gam*Gsel*slsep*t*spos;

rfmb = sum( repmat(rfsb(:),[1 mb]).* exp(1i*phi_sel) ,2);

% Convert from rad to mT using the dwell-time that scales it bmax
dtx = max(abs(rfmb))/(gam*maxb1);
rfmb = rfmb/(gam*dtx);

%%
if halfShift
  rfSoft = rfSoft.*((-1).^[0:Npulses-1]);
end

% design the blip trapezoid. Units irrelevant for dotrap fn.
garea = 1/slsep/(gam/2/pi); % mT/m*s, k-area of each blip
gzblip = dotrap(garea,gmax,gslew,dt).'; 
if ~minRFdur
    % matched-duration RF subpulses
%     hpw = ceil(max(abs(rfSoft))./(2*pi*gambar*maxb1*dt));
    
    
    % 19/01/17 sas change the below such that it design for mixing-ratio
    % PINS pulses
    if Mixing_ratio == 1
        hpw = ceil(max(abs(rfSoft))./(gam*maxb1*dt));
    else
        hpw = ceil(max(abs(rfSoft))./(gam*maxb1/(1-Mixing_ratio)*dt));
    end

    % interleave RF subpulses with gradient blips to form full pulses
%     rf = kron(rfSoft(:),[ones(hpw,1);zeros(size(gzblip))]);

    % sas pins-rf stays the same..
    rf_pins = kron(rfSoft(:),[ones(hpw,1);zeros(size(gzblip))]);
    g_pins = repmat([zeros(hpw,1);gzblip(:)],[Npulses 1]);
    
    % sas MB is played during gradient blips
    N = ceil(length(rfmb)*dtx/dt); 
    rfmb = interp1(linspace(0,1,length(rfmb)),rfmb,linspace(0,1,N));

%   sas 17/01/17 Think this verse can do it.
    % Since rfmb is interpolated, it's BW has changed. So need to
    % re-evaluate for Gselection
    
    Gsel2 = 2*pi*tb/(N*dt)/(gam*slthick);
    rfmb_v = verse_arbg(rfmb,Gsel2*ones(N,1),dt,g_pins,dt,gslew);
    nsn = find(isnan(rfmb_v));
    if any(nsn)
        fprintf('Found isnan at in MultiPINS RF at %d elements\n',length(nsn));
        rfmb_v(nsn)=0;
    end

else

    % matched-amplitude RF subpulses
    rf_pins= [];g_pins = [];
    for ii = 1:Npulses

        % Find non-integer (ni) number of subpulse-samples
%         hpw_ni = abs(rfSoft(ii))/(gam*maxb1*dt);
%         hpw = ceil(abs(rfSoft(ii))./(gam*maxb1*dt));

        if Mixing_ratio == 1            
            hpw_ni = abs(rfSoft(ii))/(gam*maxb1*dt);
            hpw = ceil(abs(rfSoft(ii))./(gam*maxb1*dt));
        else            
            hpw_ni = abs(rfSoft(ii))/(gam*maxb1*(1-Mixing_ratio)*dt);
            hpw = ceil(abs(rfSoft(ii))./(gam*maxb1/(1-Mixing_ratio)*dt));
        end
        
        beta =maxb1/hpw * (hpw-hpw_ni);
        rf_pins = [rf_pins; sign(rfSoft(ii))*gam*dt*(maxb1-beta)*ones(hpw,1);zeros(size(gzblip))];
        g_pins = [g_pins; zeros(hpw,1); gzblip(:)];
    end   
    
    
        % sas MB is played during gradient blips
    N = ceil(length(rfmb)*dtx/dt); %<--
    rfmb = interp1(linspace(0,1,length(rfmb)),rfmb,linspace(0,1,N));

%   sas 17/01/17 Think this verse can do it.
    % Since rfmb is interpolated, it's BW has changed. So need to
    % re-evaluate for Gselection
    
    Gsel2 = 2*pi*tb/(N*dt)/(gam*slthick);
    rfmb_v = verse_arbg(rfmb,Gsel2*ones(N,1),dt,g_pins,dt,gslew);
    nsn = find(isnan(rfmb_v));
    if any(nsn)
        fprintf('Found isnan at in MultiPINS RF at %d elements\n',length(nsn));
        rfmb_v(nsn)=0;
    end
        
end
% rf = rf./(sum(rf(:))*2*pi*gambar*dt)*sum(rfSoft); % convert to gauss
rf_pins = rf_pins./(sum(rf_pins(:))*gam*dt)*sum(rfSoft); % convert to mT

rf = (1-Mixing_ratio)*rf_pins + (Mixing_ratio*rfmb_v);
g = g_pins;

% remove the last blip, if 'se'
g = g(1:end-length(gzblip)); % flip sign (gradient reversal to suppress fat) and remove last blip
rf = rf(1:end-length(gzblip));
