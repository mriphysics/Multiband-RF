function [ B1_demand, B1_actual,G_demand,G_actual,B1T_store] = dz_MBverseGirf(rf_init,Gz,dt,maxg,maxgslew,b1max,verse_singleband,mb,Gsel_gap,dt_os,epsilon,max_iterations,AM_only,girf)
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
% Gsel_gap : Selection gradient [mT/m] * Gap [m]. Need this for evaluating
%           modulation function when verse_singleband = 1.
% dt_os  : Time over-sampling factor for increasing time-resolution. 
%          Recommended: For verse_singleband = 0 set this to 1 (i.e. no oversampling) as
%          it likely causes interpolation errors. For verse_singleband = 1
%          set this to 2-5 depending on slice-separation.
% epsilon: Double [mT]. A variable to design a VERSE RF pulse for
%          b1max-epsilon. Recommended: ~2e-4 for MBv, 1e-5 for vMB.
% max_iterations  : Stop after a number of iterations.
% AM_only: Only matters for vMB pulses. Set to 0 for Complex-valued MB pulse (Wong ISMRM 2012). Set to 1 for AM-only MB pulse (Malik ISMRM 2013)
% girf   : If this is a measured GIRF, input as a structure such that the function
%           gradient_distort_GIRF can useit. Else, for a model GIRF
%           (exponential) set this to a 1x2 vector. 
%         Recommended: [42 42]*1e-6.

    if ~isvector(rf_init) || ~isvector(Gz)
        error('RF and Gradient must be vectors');
    end
        
    if nargin < 14
        girf = load('GIRF_all_axes_20140729');
        warning('Defining GIRF internally');
        measured_girf = 1;
    end
    
    if isstruct(girf)
        measured_girf = 1;
        fprintf('Using a measured GIRF\n');
    else
        measured_girf = 0;
        fprintf('Using a modelled GIRF with tau = %.2f us\n',girf(1)*1e6);
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
        
        % 06-Jul 2018 Insert case for no phase-optimization.
        if AM_only == -1
            phi_sol_PO = zeros(size(phi_sol_PO));
        end
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
    ii = 0;
    current_max = 1;
    while and(current_max>b1max,ii<max_iterations)
        ii = ii +1;

        if ii ==1 
% % initialise Alpha
            alpha = ones(1,length(G0));  
            dt0 = dt;
            t0 = t;          
        else                        
            if verse_singleband
                verse_in.b = rfsb_actual(:)*10; % Convert from mT --> Gauss
            else
                verse_in.b = B1_actual(:)*10; %<--- Multiply with alpha internally. Convert from mT --> Gauss
            end
%             verse_in.bmax = (b1max-epsilon)*10;
%             if ii>50
%                 verse_in.bmax = (b1max-ii/10*epsilon)*10;
%             end

% Might try having lower b1 from the 2nd iteration onwards.
            verse_in.bmax = (verse_in.bmax - epsilon*10);

            alpha = ones(Ntv,1); %<--- now a redundant variable. Might remove in later version.
            Gn = VERSE_out.g; %<--- Already in Gauss, no conversion needed.
            verse_in.Gmod = sum((Gn).^2,2).^0.5; %<--- Already in Gauss
            dt = dtv;
            t = dt*(0:Ntv-1)';
            k0 = gamma_mT*cumtrapz(t,Gn*10);%<-- Gauss/m to mT/m
            Kcm = k0 / (2*pi*100);% back to 1/cm 
            
        end

        % Run time-optimal VERSE
        [C,time,g,slew,k, phi, sta, stb,p_of_t,VERSE_out] = minTimeGradient_VERSE_girf(Kcm,0,0,Gmax,Smax,dt0*1e3/dt_os,[],0,verse_in,alpha);

        B1_demand = VERSE_out.bt(:)/10; %<-- Gauss -> mT
        G_demand = VERSE_out.g*10;  %<-- Gauss -> mT/m
        dtv = 1e-3*(VERSE_out.t(2)-VERSE_out.t(1));
        Ntv = length(B1_demand);      
        
        % Include girf to obtain G_actual
%         % Define GIRF externally
%         girf = load('GIRF_all_axes_20140729');
        if measured_girf
            G_actual = gradient_distort_GIRF(G_demand,girf.ff,girf.Hw,dtv,500);
        else %<-- Use mono-exponential GIRF
%             G_actual = gradient_distort(G_demand,girf(1),girf(2)*ones(1,3),dtv,500);
            G_actual = gradient_distort_FT(G_demand,girf(1),girf(2)*ones(1,3),dtv,500);
        end

        if verse_singleband
        % Evaluate rfsb_actual
            rfsb_actual = verse_arbg(rf_init,G0(:,3),dt0,G_actual(:,3),dtv);
            if any(isnan(rfsb_actual))
                keyboard;
            end
		% To evaluate the modulation for altered excitation k-space, first
		% design the modulation function for a linear excitation k-space
		% and reshape it based on the verse gradient.
        
            
    %%        
            
            spos = (1:mb)-(mb+1)/2; 
            phi_sel = gamma_mT*Gsel_gap*t0(:)*spos;

		% Not summed
            mod_normal = exp(1i*phi_sel + repmat(1i*phi_sol_PO,[nt 1]));    

%         %     % Test MB of normal MB pulse.
%             FOV = 0.12;
%             Nz = 2048;
%             z = linspace(-FOV/2,FOV/2,Nz)';
%             pos = [z(:)*0 z(:)*0 z(:)];
%             pulseMB_normal = sum( repmat(rf_init,[1 mb]).*mod_normal,2);
%             [~,~,~,~,an,bn] = blochsim_CK(pulseMB_normal,G0,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dt);
% %             mxyn = 2*conj(an(:,end)).*bn(:,end);
%             mxy0 = bn(:,end).^2;
%             figure;
%             
%             plot(z,abs(mxy0));
%             dum = 1;
            
%%
        % Correct each modulation function separately.      
            clear mod_actual mod_demand %<-- clear versed modulation function from prev iteration.            
            for i = 1:mb
                mod_actual(:,i) = verse_arbg(mod_normal(:,i),G0(:,3),dt0,G_actual(:,3),dtv);
                mod_demand(:,i) = verse_arbg(mod_normal(:,i),G0(:,3),dt0,G_demand(:,3),dtv);
            end
            dum = 1;
            
         % SAS 12/01/17 Scale modulation function such that each component has unit
         % amplitude. See notes for illustration.
            mod_actual = mod_actual./abs(mod_actual);
            mod_actual(isnan(mod_actual))=0;
            mod_demand = mod_demand ./abs(mod_demand);
            mod_demand(isnan(mod_demand))=0;

            B1_actual = sum( repmat(rfsb_actual,[1 mb]).*mod_actual,2);
            
            % 05/09/17 Should be the same as adding MB slice-shifted versions
            % of the single-band version, underying a ktx modulation
            % function. However not possible in this function because
            % Gsel_gap [mT] cannot be converted back into gap [m].
            
%             tvmb = (0:(Ntv-1))*dtv;
%             ktx = gamma_mT*flipud(cumtrapz(tvmb,flipud(G_actual(:,3))));
%             for i = 1:mb
%                 mod_actual_ktx(:,i) =exp(1i*ktx*0.028*spos(i)+repmat(1i*phi_sol_PO(i),[Ntv 1]));
%             end
            % Debug for B1_actual                        
%             [~,~,~,~,an,bn] = blochsim_CK(B1_actual,G_actual,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dtv);
%             mxyn = 2*conj(an(:,end)).*bn(:,end);
%             mxyn = bn(:,end).^2;
%             figure;
%             hold on;            
%             plot(z,abs(mxyn));
%             dum = 1;
            
            % Debug for B1_actual_ktx 
%             B1_actual_ktx =sum( repmat(rfsb_actual,[1 mb]).*mod_actual_ktx,2);             
%             [~,~,~,~,an,bn] = blochsim_CK(B1_actual_ktx,G_actual,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dtv);
%             mxyn_ktx = bn(:,end).^2;
%             plot(z,abs(mxyn_ktx));
%             dum = 1;                       
            
%             B1_demand_MB = sum( repmat(B1_demand,[1 mb]).*mod_actual,2);            
            % Note that B1_demand in this method is a VERSE'd Single-band
            % waveform. Hence a distinction exists between B1_demand_MB and
            % B1_demand. But I'm not sure why it's worth keeping a
            % Singleband version of B1_demand so just replace it.
            
            % 13/02/17 - Changed my mind on this. What's interesting to
            % return is B1_demand made from the "demand" modulation
            % function, which does not make use of any GIRF.
            B1_demand_MB = sum( repmat(B1_demand,[1 mb]).*mod_demand,2);
            B1_demand = B1_demand_MB; clear B1_demand_MB;

%             [~,~,~,~,an,bn] = blochsim_CK(B1_demand,G_demand,pos,ones([Nz 1]),zeros([Nz 1]),'dt',dtv);
%             mxyn = 2*conj(an(:,end)).*bn(:,end);
%             figure;
%             hold on;            
%             plot(z,abs(mxyn));
%             legend('MB const','B1_{act} G_{act}','B1_{dem} G_{dem}');
%             dum = 1;  
            rfsb_actual_store{ii} = rfsb_actual;
        else
            % Evaluate RF_actual
            B1_actual = verse_arbg(rf_init(:),G0(:,3),dt0,G_actual(:,3),dtv);            
        end
        
        current_max = max(abs(B1_actual));
        
        dtv_store{ii} = dtv;
        G_demand_store{ii} = G_demand;
        G_actual_store{ii} = G_actual;
        B1_demand_store{ii} = B1_demand; 
        B1_actual_store{ii} = B1_actual;
        B1T_store{ii} = current_max*dtv*Ntv*1e6;
        fprintf('it nr:%.d /%.d B1:%.4f B1*T:%.4f Nt: %.d Ntv:%.d \n',ii,max_iterations,current_max*1e3,current_max*dtv*Ntv*1e6,nt,Ntv)
        
%  % Debug
%     fh = figure;
%     subplot(311);plot(t0,abs(pulseMB));hold on;plot(t0,G0(:,3)/1e3);
%     title('Const G');
%     subplot(312);plot(tv,abs(B1_demand));hold on;plot(tv,G_demand(:,3)/1e3);
%     title('B1_{dem}, G_{dem}');
%     subplot(313);plot(tv,abs(B1_actual));hold on;plot(tv,G_actual(:,3)/1e3);
%     title('B1_{act}, G_{act}');
%     set(fh,'pos',[   680   139   656   839]);
%     dum = 1;
    end
%     if verse_singleband
    if 0
        dum = 1;

%         %%
        its = 1:ii;
        for l = its
%             RF = rfsb_actual_store{l};
            RF = B1_actual_store{l};
            teff(l) = length(RF)*dtv_store{l}*max(abs(RF));            
        end
        figure;plot(its,teff);
        %%
%         for l = [3 5 10]
        for l = [1 2 3]
            fh = figure;
            set(fh,'pos',[[680 706 1228 272]]);
            dtv = dtv_store{l};
            Nv = length(rfsb_actual_store{l});
            t = 0:dtv:(Nv-1)*dtv;
            subplot(131);plot(t,rfsb_actual_store{l}*1);            
            subplot(132);plot(t,G_demand_store{l}(:,3));title(sprintf('G_{dem},it:%d',l));
            hold on;plot(t,G_actual_store{l}(:,3));
%             legend('dem','act');
            subplot(133);
            plot(t,G_demand_store{l}(:,3)-G_actual_store{l}(:,3));
            title('G_{dem} - G_{act}');
            axis([-inf inf -0.2 0.2]);
        end
end

