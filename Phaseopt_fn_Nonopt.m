function [ rfmb ] = Phaseopt_fn_Nonopt(rfsb,mb,tb,bs)
%PHASEOPT_FN Summary of this function goes here
%   Detailed explanation goes here
% 24/11/16 - Function which designs non-optimized MB waveform, for testing.

N = length(rfsb);
t = 0:1/(N-1):1;

spos = (1:mb)-(mb+1)/2; 
phi_sel = 2*pi*tb*bs*t(:)*spos;

phi_sol_PO = zeros(1,mb); 

rfmb = sum( repmat(rfsb(:),[1 mb]).* exp(1i*phi_sel + repmat(1i*phi_sol_PO,[N 1])) ,2);

end

