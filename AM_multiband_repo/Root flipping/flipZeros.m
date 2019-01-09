function flip = flipZeros(N,r,wp,bp)

% Beta polynomial root flipping 
% Inputs:
%   b: beta coefficients
%   ntrials: # monte-carlo trials
%   bc: normalized band centers in radians (-pi to pi)
%   flip: flip angle in radians
%   d1: passband ripple
%   TBW: time-bandwidth product of pulse
% Outputs:
%   rfout: rf with lowest peak amplitude
%   bestroots: roots corresponding to rfout
%   bestflips: flipped root pattern
%   initroots: initial roots
% 
% Adapted from Miki Lustig's flipZeros code by Will Grissom and Anuj Sharma
% Vanderbilt University, 2015

% find the indices of the bands
w = (-N/2+1:N/2)/N*2*pi;
idxPass = [];
for ii = 1:length(bp)
    idxPass = [idxPass find(w >= (bp(ii)-wp) & w <= (bp(ii)+wp))];
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

flip = xor(flip,tmp);   % apply opposite flipping pattern to other half of complex plane
                        % - if there are stopband roots not on the unit circle, they will be flipped inadvertently



