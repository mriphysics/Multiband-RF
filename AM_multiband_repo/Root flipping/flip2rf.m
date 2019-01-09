function [rf180,peak,doflip] = flip2rf(N,r,wp,bp,p,symtype,d1)

% 07/09/2015 sas - function that takes in as input a binary vector p of
% length, the number of passbands (or half if symmetric constraint is
% imposed) and returns a rf pulse.

% In my code, the rf pulse will be the refocusing pulse for which a
% excitation pulse will be generated for.

% Adapted code from Grissoms flipZeros function, which in turn was adapted
% from Miki Lustig. This function will be used in a genetic algorithm
% implementation of root-flipping

% Please use under MIT license (Copyright (c) 2015 mriphysics)
% Samy Abo Seada, Anthony Price, Jo Hajnal and Shaihan Malik. January 2017
% King's College London

w = (-N/2+1:N/2)/N*2*pi;
idxPass = [];
for ii = 1:length(bp)
    idxPass = [idxPass find(w >= (bp(ii)-wp) & w <= (bp(ii)+wp))];
end
nPb=length(idxPass);
% distribute probability across passband roots
pFlip = zeros(1,length(r));
if strcmp(symtype,'TS')==1 %|| strcmp(symtype,'CS')==1
    pFlip(idxPass) = [p,zeros(1,ceil(nPb/2))];
elseif strcmp(symtype,'CS')==1
    pFlip(idxPass) = [p,zeros(1,ceil(nPb/2))];
    
        % 26/10 Find the element related to final bp(end)+wp
    fnz = floor( (N/2 + 1) + (bp(1)-wp)*N/2/pi );
    % Find all roots on the bottom-halve of the complex plane that are beyond
    % the final flipped- passband root.
    tmp = and(angle(r)<0 , angle(r)< angle((r(fnz))) );
    % 26/10/15 Sinusoidal stop-band root-flipping
    dt = -3*pi:6*pi/(fnz-2):3*pi;
    tmp(1:fnz-1) = round(0.5+0.87*sin(dt));

    pFlip = or(pFlip,tmp');
else
    warning('Symmetry type not specified - Resorting to time-symmetric');
    pFlip(idxPass) = [p,zeros(1,ceil(nPb/2))];
end

% find associated roots on the other half of the plane; flip them too
flip=pFlip;

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

if strcmp(symtype,'TS')==1
    tmp = (angle(r)<0).';
    doflip = xor(flip,tmp);                        
elseif strcmp(symtype,'CS')==1
    doflip = flip;
end


rt = r;
rt(doflip == 1) = conj(1./rt(doflip==1));
    
    
    % get root-flipped RF
    R = poly(leja_fast(rt)); % get polynomial coefficients back from flipped roots
    R = R/max(abs(freqz(R))); % normalized by max of filter response
    bt = R*sin(pi/2 + atan(d1*2)/2); % scale to target flip
    
    [~,at]=b2amp(bt);

    rft=islr(at,bt);
       

% sas - call min peak refocusing pulse rf180. Assign outside the MC
% iterations.
if strcmp(symtype,'CS')
    rf180=imag(rft);
else
    rf180=rft;
end

peak=max(abs(rf180));