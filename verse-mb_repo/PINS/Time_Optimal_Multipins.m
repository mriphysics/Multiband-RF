function [ rfout,gout,Mixing_ratio,RF_energy_vs_mixingRatio] = Time_Optimal_Multipins( rfsb,mb,tb,slsep,slthick,de,maxb1,gmax,gslew,dt,minRFdur,halfShift,Mixing_ratio_vector)
% This function is a shell function around the dz_MultiPINS function to
% design Time-optimal MultiPINS pulses. The different mixing-ratio tried
% should be input the Mixing_ratio_vector.

% 23/04/2018 sas - version to release as part of verse-mb publication. 
N = length(rfsb);
k = length(Mixing_ratio_vector);
viol = zeros(k,1);
rfe_mx = zeros(k,1);
for i = 1:k
    Mixing_ratio = Mixing_ratio_vector(i);
    [rf,g] = dz_Multipins(rfsb,mb,tb,slsep,slthick,de,maxb1,gmax,...
        gslew,dt,minRFdur,halfShift,Mixing_ratio);  
    nsn = find(isnan(rf));
    if any(nsn)
%         fprintf('Found isnan at in MultiPINS RF at %d elements\n',length(nsn));
        rf(nsn)=0;
    end
    
    % Record designs where maxB1 is violated by MB component.
    if max(abs(rf))>maxb1 && i > 1
        viol(i) = 1;
%         break
    end
    
    rfe = sum(abs(rf).^2);
    % Record
    rf_mx{i} = rf;
    g_mx{i}  = g;
    rfe_mx(i) = rfe;    
    T_mx(i) = length(rf)*dt;
    
end

% Output the RF waveform which does not violate b1max with minimum duration
l = find(viol,1) -1;
if isempty(l) %<-- if none of the designs violated the b1-constraint, return the design with minimum duration.
    [~,l] = min(T_mx);
    rfout = rf_mx{l};
    gout = g_mx{l};
    Mixing_ratio = Mixing_ratio_vector(l);
else
    rfout = rf_mx{l};
    gout = g_mx{l};
    Mixing_ratio = Mixing_ratio_vector(l);
end
RF_energy_vs_mixingRatio = rfe_mx;
