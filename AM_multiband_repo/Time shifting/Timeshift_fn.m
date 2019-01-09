function [ rfmb ] = Timeshift_fn( rfsb,mb,tb,bs,frac_shift,AM_only )
%TIMESHIFT_FN - Creates a phase-optimized Time-shifted Multiband pulse
%based on a specified time-shift.

%   rfsb : Nx1 or 1x N vector specifing the underlying single-band waveform.
%   mb   : Multiband factor. The number of slices required. Integer valued.
%   tb   : Time-bandwidth prodcut [dimensionless].
%   bs   : band-separation in unit of slice-thicknesses, slice centre to
%   frac_shift: fractional shift between 0 and 1. A shift of 0 implies
%               no-shift and results in Wong's method. A shift of 1 doubles the duration.

%   AM_only: Boolean value. Set to 0 for Wong [2012] method and to 1 to
%   create a phase-optimized Multiband pulse which only contains Amplitude
%   modulation.

% Please use under MIT license (Copyright (c) 2015 mriphysics)
% Samy Abo Seada, Anthony Price, Jo Hajnal and Shaihan Malik. January 2017
% King's College London

        rfsb = rfsb(:);
        rfsb = rfsb';
        % Create a zero-padded version of the single-band waveform
        rfsb_zp = [rfsb 0*rfsb]';    

        N = length(rfsb);
        M = length(rfsb_zp);
        t = 0:1/(M-1):1;
        
        spos = (1:mb)-(mb+1)/2;
        phi_sel = 2*2*pi*tb*bs*t(:)*spos;

        if AM_only
            %%% AM only means slices are paired up
            Nvar=ceil(mb/2);
        else
            Nvar = mb;
        end

        %%% Define lower and upper-bound for constrained optimization
        lb = zeros([1 Nvar]);
        ub = [2*pi*ones([1 Nvar])];

        tau = linspace(0,frac_shift*length(rfsb),Nvar);

        %%% Define cost function here for these shifts
        if AM_only
            pfun = @(x)(timeshiftmb(rfsb_zp,phi_sel,x,tau,'AM'));
            maxfun =  @(x)(max(abs(pfun(x))));
            costfun = @(x)(max(abs(real(pfun(x))))+100*max(abs(imag(pfun(x)))));
        else
            pfun = @(x)(timeshiftmb(rfsb_zp,phi_sel,x,tau));
            maxfun = @(x)(max(abs(pfun(x))));
            costfun = maxfun;
        end
    %         dur = @(x)((max(tau)-min(tau))/(M/2));

        %%% GENETIC ALGORITHM
        options = gaoptimset;
        options = gaoptimset(options,'Display', 'off','Tolfun',1e-3);
        options = gaoptimset(options,'StallGenLimit',50);
    %                     options = gaoptimset(options,'PlotFcns', {  @gaplotbestf @gaplotbestindiv @gaplotscorediversity @gaplotscores });
        options = gaoptimset(options,'UseParallel','always');%<- use this on Arnold
        options = gaoptimset(options,'PopulationSize',50);
        % Seed one solution as the non-shifted optimum
        if AM_only
            ww = load('bmax_conj');
        else
            ww = load('bmax_wong');
        end
        xw = ww.pstore{mb}(1:Nvar);
        % Make them wrap to be within bounds (0 to 2pi)
        xw = angle(exp(1i*xw))+pi;

        % If this is not the first shift value then also add the last
        % solution to the initial population
        options = gaoptimset(options,'InitialPopulation',xw);

        [x_sol,fval,exitflag,output,population,score] = ...
            ga(@(x)(costfun(x)),Nvar,[],[],[],[],lb,ub,[],[],options);

        %%% Add in fminsearch to help get to the local minimum
        opt = optimset('MaxIter',100);
        x_sol = fminsearch(costfun,x_sol,options);

        if AM_only
            [pulseMB,~] = timeshiftmb(rfsb_zp,phi_sel,x_sol,tau,'AM');
        else
            [pulseMB,~] = timeshiftmb(rfsb_zp,phi_sel,x_sol,tau);                    
        end
        M_mb = ceil((1 + frac_shift) * length(rfsb)); %<-- This is how long pulseMB should have non-zero values up to.
        %
        rfmb = pulseMB(1:M_mb);

end

function [pulseMB,pulses_all] = timeshiftmb(rfsb,phi_sel,phi_off,tau,varargin)

M = length(rfsb);
mb = size(phi_sel,2);

%%% remove offsets from timeshifts
tau = tau - min(tau);

%%% If addition argument, then assume we have AM version which means paired
%%% slices
if ~isempty(varargin)

    %%% Mapping for indices
    if mod(mb,2)==0 
        tau = [tau fliplr(tau)];             % Each pair has SAME shift   
        phi_off = [phi_off -fliplr(phi_off)];% Phase is conjugated
    else
        tau = [tau fliplr(tau(1:end-1))];
        phi_off = [phi_off -fliplr(phi_off(1:end-1))];
    end

    
    
else
    if (length(phi_off)~=mb)||(length(tau)~=mb)
        fprintf(1,'Error: incorrect number of phase offsets/time shifts\n');
        return
    end

end
    
    
%%% Replicate the single band pulse mb times and include the selection
%%% gradient associated phase (phi_sel = gamma_mT*Gsel*gap*t*spos)

pulseMB = zeros([M mb]);

for ii=1:mb
    % First deal with selection gradient induced phase and offset applied
    pulseMB(:,ii) = rfsb .* exp(1i*(phi_sel(:,ii)+phi_off(ii)));
    
    % use linear resampling? Nearest doesn't require correction later..
    pulseMB(:,ii) = interp1(1:M,pulseMB(:,ii),(1:M)-tau(ii),'nearest');
end
pulseMB(isnan(pulseMB))=0;

if nargout==2
    pulses_all = pulseMB;
end
pulseMB = sum(pulseMB,2);

end