function [ an,ac] = b2amp( bc,osn )
% Generate a minimum-phase filter from a filter bc. This function is part
% of the Shinnar-Le Roux code for designing large-tip RF pulses.

% Up until now I used a function called GenAn_Minphi.m However, I've now
% realised that it's important to zeropad the time-profile of the input
% filter before taking the fft and the Hilbert transform. Over the last
% three weeks (16/09/15) I made detailed notes about how this affects the
% subsequent excitation profiles. Furthermore, Pauly's b2a.c code also
% zeropads by I believe 16*length(bc).

% Please use under MIT license (Copyright (c) 2015 mriphysics)
% Samy Abo Seada, Anthony Price, Jo Hajnal and Shaihan Malik. January 2017
% King's College London

N = length(bc);
% Find Over-sampling factor. The exp-term makes sure the zeropadded length
% is a power of two (for FFT efficiency), even if N is not. Run
% exp(ceil(log(1:10)/log(2))*log(2))if you're confused.

% OS = 32*exp(ceil(log(N)/log(2))*log(2));
% 18/03/16 Allow user input for oversampling

if nargin < 2
    osn = 32;
end
OS = osn*exp(ceil(log(N)/log(2))*log(2));
bc_zeropadded =...
    [zeros(ceil((OS-N-1)/2),1);
    bc(:);
    zeros(ceil((OS-N-1)/2),1)];

if max(abs(fft(bc_zeropadded)).^2)>1
    bc_zeropadded=bc_zeropadded./max(abs(fft(bc_zeropadded))).^2;
end

if mod( log2(length(bc_zeropadded)), 1) ~= 0
    warning('length bc_zeropadded not a power of two');
end

bn = fft(bc_zeropadded);

anmag=sqrt(1-abs(bn).^2);
an=anmag.*exp(1i*imag(hilbert(log(anmag))));

ac_zeropadded = fft(an)/OS;
ac = ac_zeropadded(1:N);

if iscolumn(bc) == 0
    ac = ac';
end

end