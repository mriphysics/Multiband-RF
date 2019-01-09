function [ rf90, rf180,tb] = rootflip_fn(n,nb,tblin,bandsep,CS)
% 31/10/2015 sas Function that returns a pair of root-flipped pulses for a
% spin-echo sequence. The returned pulses are scaled in radians and defined
% in normalized time given by n-time points for the refocusing pulse and
% 2*n-time points for the excitation pulse.

% Input:
% n  : number of time points
% nb : number of bands/slices
% tblin : Time-bandwidth product
% bandsep : number of slices of separation between bands (measured from
% passband centre to passband centre)
% alignedecho = 0; % 0 = design a minimum-duration 90
% CS :% Boolean for conjugate Symmetry. 1 for AM pulses

% Output:
% rf90 - excitation pulse in radians, with length 2*n
% rf180- refocusing pulse in radians, with length n
% tb   - the new time-bandwidth product after the lin-phase to min-phase
% transformation. see eq.10 in Sharma et al. Use this outside the function
% to find the Bandwidth and hence gradient waveform.

% Adapated from code released by Sharma et al
% (http://www.vuiis.vanderbilt.edu/~grissowa/software.html)

% Please use under MIT license (Copyright (c) 2015 mriphysics)
% Samy Abo Seada, Anthony Price, Jo Hajnal and Shaihan Malik. January 2017
% King's College London

bandsep = bandsep*tblin;

d1 = 0.01;             % combined Mxy ripple, passband
d2 = 0.01;             % combined Mxy ripple, stopband
                        % 1 Gauss = 0.1mT, or 1e-4T
rftype = 'matched';     % 'matched' or '2xrefd' (twice-refocused)

osfact = 10; % oversampling factor

N = 2*(n-1); % length of linear-phase filter we will factor to get min-phase filter

% directly design a MB min-phase filter, whose roots we can flip
% sas - adjust the ripples depending on 'matched' or 'twice-refocused'
if strcmp(rftype,'matched')
  d1 = d1/4; % target beta passband ripple, considering full 90-180 pair
  d2 = (d2/sqrt(2))^0.25; % target beta stopband ripple, considering full 90-180 pair
else
  d1 = d1/8; % target beta passband ripple, considering twice-refocused
  d2 = d2.^(1/4); % target beta stopband ripple, considering twice-refocused
end

nn = (0:N/2*osfact)'/(N/2*osfact);  % 0 to 1 - profile indices
d = zeros(N/2*osfact+1,1);          % passband mask
s = zeros(N/2*osfact+1,1);          % stopband mask
wts = zeros(N/2*osfact+1,1);        % ripple taper weights

dinfmin = 1/2*dinf(2*d1,d2^2/2); % d-infinity for a min-phase pulse with beta ripples (d1,d2)
dinflin = dinf(d1,d2);      % d-infinity for a linear phase pulse with the same ripples


tb = tblin/dinflin*dinfmin; % scale TBW product so as to get the same transition 
                            % width as linear phase pulse with same ripples, 
                            % after scaling back to desired slice thickness. This 
                            % makes comparison to other MB excitations more 
                            % meaningful, since all will have same slice characteristics.
w = dinfmin/tb; % transition width

if rem(nb,2) % if Nb odd
    % start out the f, m and w vectors with the DC band
    f = [0 (1-w)*(tb/2) (1+w)*(tb/2)];%*di/dilp;
    d = nn <= f(2)/(n/2); % target pattern
    wts = 1./abs(nn).^2; % quadratically-decaying ripple weights
else
    f = 0;
end   

% add non-DC bands to the profiles
for ii = 1:floor(nb/2)
    cent = (ii - (rem(nb,2) == 0)/2)*(bandsep)*dinfmin/dinflin;
%     Old version of code used:
%     cent = (ii - (rem(nb,2) == 0)/2)*(bandsep-2/osfact)*dinfmin/dinflin
    f = [f (cent-(1+w)*(tb/2)) (cent-(1-w)*(tb/2)) (cent+(1-w)*(tb/2)) (cent+(1+w)*(tb/2))];
    d = d | (nn >= f(end-2)/(n/2) & nn <= f(end-1)/(n/2));
    s = s | (nn >= f(end-4)/(n/2) & nn <= f(end-3)/(n/2));
    nnc = nn - (f(end-1)+f(end-2))/2/(n/2); % indices centered with passband for weight calcs
    wts = max(wts,1./abs(nnc).^2); % quadratically-decaying ripple weights
end
% append the last stopband
s = s | (nn >= f(end)/(n/2));
wts = wts.*s;wts(isnan(wts)) = 0;wts = wts./max(wts);

% build system matrix for cvx design
A = 2*cos(2*pi*(0:N/2*osfact)'*(-N/2:0)/(N*osfact));A(:,end) = 1/2*A(:,end);

% mask everything to get rid of transition bands
% Note how Ad has fewer rows than A, because the transition bands are not
% included.
Ad = A(s | d,:);
dd = double(d(s | d));
ss = wts(s | d).*double(s(s | d));

% use cvx to do the constrained optimization
cvx_begin
  variable delta(1) 
  variable x(N/2+1)
  minimize( delta )
  subject to
    -delta*dd <= Ad*x - dd <= delta*dd + delta*d2^2/(2*d1)*ss
cvx_end

% stack two halves together to get full linear-phase filter
x = [x;x(end-1:-1:1)]';
blin=x;
% factor the linear phase filter to get a min-phase filter b
% b = real(fmp(x));
% sas - Use a new function fmp2 instead of fmp to distinguish from the fmp
% function in Pauly's RF toolbox, which has lower oversampling can
% introduce amplified ripples.
b = real(fmp2(x));
b = b(end:-1:1);

% root flip h!
ntrials = ceil(tb)*nb*100; % number of monte carlo trials

% sas- Note How TBW scaling is used in the selection of passband locations
bp = 2*((1:nb)-1/2-nb/2)*bandsep/n*pi*dinfmin/dinflin; % centers of passbands

flip=pi;
N = length(b); % number of time points in pulse

b = b./max(abs(fft(b))); % normalized beta coefficients

b = b*sin(flip/2 + atan(d1*2)/2); % scale to target flip angle

[~,a] = b2amp(b);

rfinit_max = max(abs(islr(a,b))); % get initial peak RF amplitude
r = roots(b); % calculate roots of beta polynomial

% sort the roots by phase angle
[~,idx] = sort(angle(r)); 
r = r(idx);

btmx = zeros(ntrials,n);

minpeakrf = Inf;
maxpeakrf = rfinit_max;

% Try using GA tool-box. If not available on system, resort to simple
% Monte-carlo.
try  
    
%     error('Test - Default to MC optimization');
    fprintf('Optimizing Root-pattern using Genetic algorithm\n');
    nPb=nrPb_fn(N-1,(tb/N*pi)*3,bp-pi/N); % Find number of pass-band roots
    if CS==0
        symtype='TS';
    elseif CS==1
        symtype='CS';
    end
    % Number of optimization variables
    Nvars=floor(nPb/2);

    costfun=@(p)flip2peak(N-1,r,(tb/N*pi)*3,bp-pi/N,p,symtype,d1);

    options = gaoptimset;
    options = gaoptimset(options,'Display', 'off','Tolfun',1e-6);
    options = gaoptimset(options,'UseParallel','always');
    % options = gaoptimset(options,'PopulationSize',100);
    % 25 Jan dataset options:

    options = gaoptimset(options,'Display', 'off','Tolfun',1e-6,'TolCon',1e-3);

    popt=ga(costfun,Nvars,[],[],[],[],zeros(1,Nvars),ones(1,Nvars),[],1:Nvars,options);
    [rf180,optpeak,bestflips] = flip2rf(N-1,r,(tb/N*pi)*3,bp-pi/N,popt,symtype,d1);

    fprintf('Peak went from %d to %d rad\n',rfinit_max,optpeak);

    % 18/04/16 Undo the -conj() in the islr1 which happens inside flip2rf:

    N = length(rf180); % # time points
    [~,b90di] = rf2ab_sas(rf180,(-N/2:1/2:N/2-1/2)',0); % get beta of 180 pulse

    gridx=(-N/2:1/2:N/2-1/2)';

    b90d = (b90di.^2)/sqrt(2); % target 90-deg beta profile
    b90d= -conj(b90d);

    bx=fftshift( fft(ifftshift(b90d))/length(b90d) );

    [~,ax] = b2amp(bx);

    rf90 = -1i*conj(islr(ax,bx));
    if CS==1
        rf90=real(rf90);
    end

catch ME
    
    fprintf('Error in Genetic-algorithm optimization. Using Monte-carlo optimization\n');
    for ii = 1:ntrials
    
    % determine which indices to flip
    if CS==1
        doflip = flipZerosCS(N-1,r,(tb/N*pi)*3,bp-pi/N);

    elseif CS==2
        mask=0.5*ones(1,N-1);
        doflip= flipAllZeros(N-1,r,(tb/N*pi)*3,bp-pi/N,mask);
    else
        doflip = flipZeros(N-1,r,(tb/N*pi)*3,bp-pi/N);
    end
% Change to flip pattern which gives conjugate symmetry
    % flip those indices
    rt = r;
    rt(doflip == 1) = conj(1./rt(doflip==1));

    % get root-flipped RF
    R = poly(leja_fast(rt)); % get polynomial coefficients back from flipped roots
    R = R/max(abs(freqz(R))); % normalized by max of filter response
    bt = R*sin(flip/2 + atan(d1*2)/2); % scale to target flip
    btmx(ii,:) = bt;
    if max(abs(fft(bt)).^2)>1
        bt=bt./max(abs(fft(bt))).^2;
    end
    
    [~,at]=b2amp(bt);
    rft   = islr(at,bt);

%     save result if better
    if max(abs(rft)) < minpeakrf
        minpeakrf = max(abs(rft));
        rfout = rft;
    end
    
    if rem(ii,50) == 0
        fprintf('Iteration %d of %d. Peak RF: %0.2d rad. Peak init RF: %0.2d rad.\n', ...
            ii,ntrials,minpeakrf,rfinit_max);        
    end
    
    end


    if CS==1
        rf180=imag(rfout);
    else
        rf180=rfout;
    end

    N = length(rf180); % # time points

    [a90,b90di] = rf2ab_sas(rf180,(-N/2:1/2:N/2-1/2)',0); % get beta of 180 pulse

    b90d = (b90di.^2)/sqrt(2); % target 90-deg beta profile
    b90d = -conj(b90d);

    bx=fftshift( fft(ifftshift(b90d))/length(b90d) );


    [~,ax] = b2amp(bx);

    rf90 = -1i*conj(islr(ax,bx));
    if CS==1
        rf90=real(rf90);
    end
    

end

% Output columns
rf90 = rf90(:);
rf180 = rf180(:);