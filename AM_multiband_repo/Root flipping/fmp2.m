%
%  Generate an equal ripple minimum phase filter starting with a linear
%  phase filter.
%  
%  Called by dzmp.

%  written by John Pauly, 1992
%  (c) Board of Trustees, Leland Stanford Junior University

% 2014 Modified by Will Grissom to use a higher oversampling factor for root-flippled multiband refocusing pulse design

% 21/9/2015 Samy Abo Seada - renamed to fmp2 such that it doesn't conflict with
% Pauly's fmp function

function hmp = fmp2(h)

l = length(h);
if rem(l,2) == 0,
   disp('filter length must be odd');
   return;
end;
lp = 128*exp(ceil(log(l)/log(2))*log(2));
hp = [zeros(1,ceil((lp-l)/2)) h zeros(1,floor((lp-l)/2))];
hpf = fftc(hp);
%hpfs = hpf-min(real(hpf));
hpfs = hpf-min(real(hpf))*1.000001;
hpfmp = mag2mp(sqrt(abs(hpfs)));
hpmp = ifft(fftshift(conj(hpfmp)));
hmp = hpmp(1:(l+1)/2);
%keyboard

