function flip = flipZerosCS(N,r,wp,bp)

% Beta polynomial root flipping for Conjugate Symmetry (needed for AM-only
% pulses).

% Inputs:
%   b: beta coefficients
%   ntrials: # monte-carlo trials
%   bc: normalized band centers in radians (-pi to pi)
%   flip: flip angle in radians
%   d1: passband ripple
%   TBW: time-bandwidth product of pulse
% Outputs:

 
% Adapted from Miki Lustig's flipZeros code by Will Grissom and Anuj Sharma
% Vanderbilt University, 2015

% Please use under MIT license (Copyright (c) 2015 mriphysics)
% Samy Abo Seada, Anthony Price, Jo Hajnal and Shaihan Malik. January 2017
% King's College London

% find the indices of the bands
w = (-N/2+1:N/2)/N*2*pi;
idxPass = [];
% idxStop = []; % sas 
for ii = 1:length(bp)
    idxPass = [idxPass find(w >= (bp(ii)-wp) & w <= (bp(ii)+wp))];
%     idxStop = [idxStop find(w < (bp(ii)-wp) || w > (bp(ii)+wp))];
end

% distribute probability across passband roots
p = 1/length(idxPass(:)):1/length(idxPass(:)):1;
pFlip = zeros(1,length(r));
pFlip(idxPass) = p;

% get flips for half of complex plane; other half will have conjugate
% flip pattern

tmp = (angle(r)<0).';
% find associated roots on the other half of the plane; flip them too
flip = [rand(1,floor(N/2))<2*pFlip(1:floor(N/2)),zeros(1,ceil(N/2))];

% 26/10/15 find the element related to final bp(end)+wp
fnz = floor( (N/2 + 1) + (bp(1)-wp)*N/2/pi );
% Find all roots on the bottom-halve of the complex plane that are beyond
% the final flipped- passband root. This is to try and aim each stop-band
tmp = and(angle(r)<0 , angle(r)< angle((r(fnz))) );

% 26/10/15 Sinusoidal stop-band root-flipping!
dt = -3*pi:6*pi/(fnz-2):3*pi;
tmp(1:fnz-1) = round(0.5+0.87*sin(dt));


% logical OR with flip
flip = or(flip,tmp');
for n=1:floor(N/2)
    if flip(n)==0
        continue;
    end
    if mod(N,2)
        idx = find(conj(r(n))==r);
    else
        find(abs(conj(r(n))-r)<1e-10);
    end
    if length(idx)>1
        disp('warning decrease epsilon');
    end
    if length(idx)<1
        disp('single root');
    end
    if length(idx)==1
        flip(idx) = 1;
    end
end
