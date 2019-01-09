function [ B1_demand,G_demand,G_actual] = dz_MBverse(rf_init,Gz,dt,maxg,maxgslew,b1max,verse_singleband,mb,slicegap,dt_os,AM_only,girf)
% Designs a Multiband VERSE RF pulse with time-optimal VERSE and
% incorporate a Gradient Impulse Response Function.

% rf_init: Nt x 1 [mT] Initial single-band or Multiband RF pulse
% Gz     : Nt x 1 [mT/m] Initial Gradient waveform associated with rf_init
% dt     : double [s] dwell-time associated with rf_init and Gz.
% maxg   : double [mT/m] Maximum Gradient amplitude to design for. ~40mT/m
% maxgslew: double [mT/m/s] Maximum gradient slew-rate. ~2e5 mT/m/s
% b1max  : double [mT] maximum peak-ampltiude of B1-waveform.
% verse_singleband : logical [0,1]. 
%           Set to 1 to design a vMB pulse. rf_init should be a single-band
%           pulse. This will be VERSE'd and then modulated to frequency
%           location according to the VERSE gradient. 

%           Set to 0 to design a MBv pulse. rf_init should be a multiband pulse. This will then be
%           versed, corrected for GIRF and reiterated until the
%           B1-overshoots are below peak B1.

% mb     : Multiband factor
% slicegap : * Gap of slices, centre-centere [m]. Need this for evaluating modulation function when verse_singleband = 1.
% dt_os  : Time over-sampling factor for increasing time-resolution. 
%          Recommended: For verse_singleband = 0 set this to 1 (i.e. no oversampling) as
%          it likely causes interpolation errors. For verse_singleband = 1
%          set this to 2-5 depending on slice-separation.
% epsilon: Double [mT]. A variable to design a VERSE RF pulse for
%          b1max-epsilon. Recommended: ~2e-4 for MBv, 1e-5 for vMB.
% AM_only: Only matters for vMB pulses. Set to 0 for Complex-valued MB pulse (Wong ISMRM 2012). Set to 1 for AM-only MB pulse (Malik ISMRM 2013)
% girf   : If this is a measured GIRF, input as a structure such that the function
%           gradient_distort_GIRF can useit. Else, for a model GIRF
%           (exponential) set this to a 1x2 vector. 
%         Recommended: [42 42]*1e-6.

if ~isvector(rf_init) || ~isvector(Gz)
    error('RF and Gradient must be vectors');
end

if isstruct(girf)
    measured_girf = 1;
    fprintf('Using a structural GIRF\n');
else
    measured_girf = 0;
    fprintf('Using an analytical GIRF with tau = %.2f us\n',girf(1)*1e6);
end
rf_init = rf_init(:);
Gz = Gz(:);
gamma_mT = 2*pi*4.257*1e4; %<--- same as in minTime gradient function

nt = length(rf_init);
t = dt*(0:nt-1)';
G0 = [0*Gz 0*Gz Gz];

verse_in = struct;
if verse_singleband

    % Load in phase-offsets
    if AM_only
        load('bmax_conj.mat')
    else
        load('bmax_wong.mat')
    end
    phi_sol_PO = cell2mat(pstore(mb));
    phi_sol_PO = angle(exp(1i*phi_sol_PO))+pi;
    % sas 12/01/17 - Initial b1max design is the required b1max divided by the predicted
    % increase in B1-amplitude after MB-modulation. This amplitude increase is
    % stored in the vector bmax.
    verse_in.bmax = ( b1max/bmax(mb) )*10;
else
    verse_in.bmax = b1max*10;
end
k0 = gamma_mT*cumtrapz(t,G0); %<--- rx k-space...

%%% Translate everything to CGS units for Lustig code
Kcm = k0 / (2*pi*100);% k-space in 1/cm
Gmax = maxg / 10; % G/cm
Smax = maxgslew / 10 /1000; % G/cm/ms

%%% RF PULSE needs to go in

verse_in.b = rf_init*10;%<-- mT->Gauss
verse_in.Gmod = sum((G0/10).^2,2).^0.5;

verse_in.os = 10; %<--- oversample factor for calculation


% Run time-optimal VERSE
[C,time,g,slew,k, phi, sta, stb,p_of_t,VERSE_out] = minTimeGradient_VERSE_girf(Kcm,0,0,Gmax,Smax,dt*1e3/dt_os,[],0,verse_in);

B1_demand = VERSE_out.bt(:)/10; %<-- Gauss -> mT
G_demand = VERSE_out.g*10;  %<-- Gauss -> mT/m
tv = 1e-3*VERSE_out.t;
dtv = 1e-3*(VERSE_out.t(2)-VERSE_out.t(1));

if verse_singleband
% To evaluate the modulation for altered excitation k-space, first
% design the modulation function for a linear excitation k-space
% and reshape it based on the verse gradient.
% 
%     spos = (1:mb)-(mb+1)/2; 
%     phi_sel = gamma_mT*slicegap*t(:)*spos;    
%     mod_normal = exp(1i*phi_sel + repmat(1i*phi_sol_PO,[nt 1]));
%     for ll = 1:mb
%         mod_demand(:,ll) = verse_arbg(mod_normal(:,ll),G0(:,3),dt,G_demand(:,3),dtv);
%     end
%  % Scale modulation function such that each component has unit amplitude. 
% 
%     mod_demand = mod_demand ./abs(mod_demand);
%     mod_demand(isnan(mod_demand))=0;
%     B1_demand_MB = sum( repmat(B1_demand,[1 mb]).*mod_demand,2);
%     B1_demand = B1_demand_MB; clear B1_demand_MB;            
%     
    %% Method 2 - evaluate with bs and slthick, using ktx.
    spos = (1:mb)-(mb+1)/2; 
    ktx = gamma_mT*flipud(cumtrapz(tv,flipud(G_demand(:,3))));
    for ll = 1:mb
        mod_demand_ktx(:,ll) =exp(1i*ktx*slicegap*spos(ll)+repmat(1i*phi_sol_PO(ll),[length(G_demand) 1]));
    end    
    B1_demand_MB =sum( repmat(B1_demand,[1 mb]).*mod_demand_ktx,2);             
    B1_demand = B1_demand_MB; clear B1_demand_MB;            
end

if measured_girf
    G_actual = gradient_distort_GIRF(G_demand,girf.ff,girf.Hw,dtv,500);
else %<-- Use mono-exponential GIRF
    G_actual = gradient_distort_FT(G_demand,girf(1),girf(2)*ones(1,3),dtv,500);
end

end


