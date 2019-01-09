function [rfmb,Gs,dtmb,tb] = multiband_rf(type,rfsb,mb,tb,bs,slthick,maxb1,maxg,maxgslew,AM_only,girf,return_gdem,gradientslopes)
%MULTIBAND_RF Summary of this function goes here
%   Detailed explanation goes here

gamma_mT = 2*pi*4.257*1e4; %<--- same as in minTime gradient function
tb =tb;
switch type
    
    case 'no' %<-- Design non-optimized multiband RF pulse.
        rfmb = Phaseopt_fn_Nonopt(rfsb,mb,tb,bs);
        dtmb=max(abs(rfmb))./(gamma_mT*maxb1);
        BW = tb/(length(rfmb)*dtmb);
        Gs = 2*pi*BW/(gamma_mT*slthick)*ones(length(rfmb),1);
        
        % Output rfmb in mT
        rfmb = rfmb/gamma_mT/dtmb;
        
        % Slope gradients
        if gradientslopes == 1
            [Gs,rfmb] = SlopeGradientZeropadRF2(Gs,rfmb,maxgslew,dtmb);
        end
    case 'po'
        rfmb = Phaseopt_fn(rfsb,mb,tb,bs,AM_only);
        dtmb=max(abs(rfmb))./(gamma_mT*maxb1);
        BW = tb/(length(rfmb)*dtmb);
        Gs = 2*pi*BW/(gamma_mT*slthick)*ones(length(rfmb),1);
        
        % Output rfmb in mT
        rfmb = rfmb/gamma_mT/dtmb;
           
        % Slope gradients
        if gradientslopes == 1
            [Gs,rfmb] = SlopeGradientZeropadRF2(Gs,rfmb,maxgslew,dtmb);
        end
    case 'ts'
        frac_shift = 0.5;
        rfmb= Timeshift_fn(rfsb,mb,tb,bs,frac_shift,AM_only); %<-- in radians.

        % Find dwell-time such to scale to b1max.
        dtmb = max(abs(rfmb))/gamma_mT/maxb1;

        % Note that the BW for a time-shifted pulse is equivalent to its
        % single-band waveform - see fig.1 in Auerbach et al. MRM 2013
        BW = tb/(length(rfsb)*dtmb);
        Gsel = 2*pi*BW/(gamma_mT*slthick);
        Gs = Gsel*ones(length(rfmb),1);
        
        % Output rfmb in mT
        rfmb = rfmb/gamma_mT/dtmb;
        
                % Slope gradients
        if gradientslopes == 1
            [Gs,rfmb] = SlopeGradientZeropadRF2(Gs,rfmb,maxgslew,dtmb);
        end
        
    case 'rf'
        [ rf90, rfmb,tb] = rootflip_fn(length(rfsb),mb,tb,bs,AM_only,'single-se');        
        dtmb = max(abs(rfmb))/(gamma_mT*maxb1);
        rfmb = rfmb(:)/(gamma_mT*dtmb);
        Gs = 2*pi*tb/(length(rfmb)*dtmb*gamma_mT*slthick) * ones(length(rfmb),1);
        
        % Slope gradients
        if gradientslopes == 1
            [Gs,rfmb] = SlopeGradientZeropadRF2(Gs,rfmb,maxgslew,dtmb);
        end
        
    case 'mbv'        
        % First design an MB pulse
        rfmb = Phaseopt_fn(rfsb,mb,tb,bs,AM_only);
        dt=max(abs(rfmb))./(gamma_mT*maxb1);
        BW = tb/(length(rfmb)*dt);
        Gs = 2*pi*BW/(gamma_mT*slthick)*ones(length(rfmb),1);
        
        % Now VERSE it.
        fprintf('Designing MB verse pulse...\n')
        verse_singleband = 0;
        dt_os = 2;
        epsilon = 1e-3;
        max_iterations = 1;

        [rfmb,~,Gs,~]= dz_MBverseGirf(rfmb/gamma_mT/dt,Gs,dt,maxg,...
        maxgslew,maxb1,verse_singleband,mb,bs*tb,dt_os,epsilon,max_iterations,AM_only,girf);
        dtv = dt/dt_os;
        % Only keep Gz
        Gs = Gs(:,3);
        dtmb = dtv;
    case 'vmb'
        % Evaluate BW of single-band pulse..
        dt = max(abs(rfsb))/(gamma_mT*maxb1);
        BW_sb = tb/(length(rfsb)*dt);
        Gsel_sb = 2*pi*BW_sb/(gamma_mT*slthick);
        Gs = Gsel_sb*ones(length(rfsb),1);
        
        % Now Apply verse MB method
        verse_singleband = 1;
        dt_os = 2;
        epsilon = 1e-4;
        max_iterations = 1;

        [rfmb, ~,Gs,~]= dz_MBverseGirf(rfsb/gamma_mT/dt,Gs,dt,maxg,...
            maxgslew,maxb1,verse_singleband,mb,Gsel_sb*bs*slthick,dt_os,epsilon,max_iterations,AM_only,girf);        
        dtv =dt/dt_os;
%         Only keep Gs
        Gs = Gs(:,3);
        dtmb = dtv;
    case 'mbvg'
        % First design an MB pulse
        rfmb = Phaseopt_fn(rfsb,mb,tb,bs,AM_only);
        dt=max(abs(rfmb))./(gamma_mT*maxb1);
        BW = tb/(length(rfmb)*dt);
        Gs = 2*pi*BW/(gamma_mT*slthick)*ones(length(rfmb),1);
        
        % Now VERSE it.
        fprintf('Designing GIRF-corrected MB verse pulse...\n')
        verse_singleband = 0;
        dt_os = 2;
        epsilon = 1e-3;
        max_iterations = 150; %<-- Should be enough.. can take a long time!
        
        [~,rfmb,Gs_demand,Gs_actual,~]= dz_MBverseGirf(rfmb/gamma_mT/dt,Gs,dt,maxg,...
            maxgslew,maxb1,verse_singleband,mb,bs*tb,dt_os,epsilon,max_iterations,AM_only,girf);

        if return_gdem
            Gs = Gs_demand;
        else
            Gs = Gs_actual;
        end
        
        dtv = dt/dt_os;
        % Only keep Gz
        Gs = Gs(:,3);
        dtmb = dtv;
    case 'vmbg'
                
        % Evaluate BW of single-band pulse..
        dt = max(abs(rfsb))/(gamma_mT*maxb1);
        BW_sb = tb/(length(rfsb)*dt);
        Gsel_sb = 2*pi*BW_sb/(gamma_mT*slthick);
        Gs = Gsel_sb*ones(length(rfsb),1);
        
        verse_singleband = 1;
        dt_os = 1;
        epsilon = 5e-5;
        
        max_iterations = 30;

        [~,rfmb,Gs_demand,Gs_actual]= dz_MBverseGirf(rfsb/gamma_mT/dt,Gs,dt,maxg,...
            maxgslew,maxb1,verse_singleband,mb,Gsel_sb*bs*slthick,dt_os,epsilon,max_iterations,AM_only,girf);        
        
        if return_gdem
            Gs = Gs_demand;
        else
            Gs = Gs_actual;
        end
        
        dtv =dt/dt_os;
%         Only keep Gs
        Gs = Gs(:,3);
        dtmb = dtv;        
                    
    case 'pins'
        dt = max(abs(rfsb))/(gamma_mT*maxb1);
        fprintf('Designing PINS pulse...\n')
        halfShift = iseven(mb); % shift pattern by 1/2 slice gap to line up with target
        mindurRF = 1;
        de = 0.01; %<-- passband ripple. 
        [rfmb,Gs] = dz_pins_arbSB(rfsb,tb,bs*slthick,slthick,de,maxb1,maxg,maxgslew,...
             dt,mindurRF,halfShift);
        dtmb = dt;
    case 'multipins'
        dt = max(abs(rfsb))/(gamma_mT*maxb1);
        halfShift = iseven(mb); % shift pattern by 1/2 slice gap to line up with target
        mindurRF = 1; % switch to use min duration RF for all subpulses
        de = 0.01; %<-- passband ripple. 
        Mixing_ratio_v = 0:0.005:1;

    [ rfmb,Gs,Mixing_ratio,~] = Time_Optimal_Multipins(...
        rfsb,mb,tb,bs*slthick,slthick,de,maxb1,maxg,maxgslew,dt,mindurRF,halfShift,Mixing_ratio_v);
        dtmb = dt;
        
end

