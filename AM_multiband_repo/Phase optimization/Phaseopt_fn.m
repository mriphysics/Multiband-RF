function [ rfmb ] = Phaseopt_fn(rfsb,mb,tb,bs,AM_only)
%PHASEOPT_FN Creates an optimized scheduled Multiband pulse based on
%pre-calculated optimal phase-offsets.

%   rfsb : Nx1 vector specifing the underlying single-band waveform.
%   mb   : Multiband factor. The number of slices required. Integer valued.
%   tb   : Time-bandwidth prodcut [dimensionless].
%   bs   : band-separation in unit of slice-thicknesses, slice centre to
%   slice centre.
%   AM_only: Boolean value. Set to 0 for Wong [2012] method and to 1 to
%   create a phase-optimized Multiband pulse which only contains Amplitude
%   modulation.

% Please use under MIT license (Copyright (c) 2015 mriphysics)
% Samy Abo Seada, Anthony Price, Jo Hajnal and Shaihan Malik. January 2017
% King's College London

N = length(rfsb);
t = 0:1/(N-1):1;

spos = (1:mb)-(mb+1)/2; 
phi_sel = 2*pi*tb*bs*t(:)*spos;

% Load in phase-offsets
if AM_only
    load('bmax_conj.mat')
else
    load('bmax_wong.mat')
end

phi_sol_PO = cell2mat(pstore(mb));
phi_sol_PO = angle(exp(1i*phi_sol_PO))+pi;

rfmb = sum( repmat(rfsb(:),[1 mb]).* exp(1i*phi_sel + repmat(1i*phi_sol_PO,[N 1])) ,2);

% If AM_only, the imaginary components of rfmb will be negligible. Better 
% to remove such that MR system does not create phase-jumps at sign changes.
if AM_only
    rfmb = real(rfmb);
end

end

