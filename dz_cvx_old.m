function [b] = dz_cvx(Nt,tb,d1,d2,quiet,phase)

    osfact = 10; % oversampling factor for frequency grid

    nn = (0:Nt/2*osfact)'/(Nt/2*osfact);  % 0 to 1 - profile indices
    d = zeros(Nt/2*osfact+1,1);          % passband mask
    s = zeros(Nt/2*osfact+1,1);          % stopband mask
    wts = zeros(Nt/2*osfact+1,1);        % ripple taper weights

%     dinflin = dinf(d1,d2);      % d-infinity for a linear phase pulse with the same ripples

%     w = dinflin/tb; % transition width
    
%     set transition width
    w = dinf(d1,d2)/tb;

%     % start out the f, m and w vectors with the DC band
%     f = [0 (1-w)*(tb/2) (1+w)*(tb/2)];%*di/dilp;
%     d = nn <= f(2)/(Nt/2); % target pattern
%     wts = 1./abs(nn).^2; % quadratically-decaying ripple weights
%     % append the last stopband
%     s = s | (nn >= f(end)/(Nt/2));
%     wts = wts.*s;wts(isnan(wts)) = 0;wts = wts./max(wts);

%-------------- debug start
    if strcmp(phase,'linear');
        % start out the f, m and w vectors with the DC band
        f = [0 (1-w)*(tb/2) (1+w)*(tb/2)];%*di/dilp;
        d = nn <= f(2)/(Nt/2); % target pattern
        wts = 1./abs(nn).^2; % quadratically-decaying ripple weights
        % append the last stopband
        s = s | (nn >= f(end)/(Nt/2));
        wts = wts.*s;wts(isnan(wts)) = 0;wts = wts./max(wts);
    else
        % start out the f, m and w vectors with the DC band
        f = [0 (1-w)*(tb/2) (1+w)*(tb/2)];%*di/dilp;
        d = nn <= f(2)/(Nt); % target pattern
        wts = 1./abs(nn).^2; % quadratically-decaying ripple weights
        % append the last stopband
        s = s | (nn >= f(end)/(Nt));
        wts = wts.*s;wts(isnan(wts)) = 0;wts = wts./max(wts);
    end
%-------------- debug end

    % build system matrix for cvx design
    A = 2*cos(2*pi*(0:Nt/2*osfact)'*(-Nt/2:0)/(Nt*osfact));
    A(:,end) = 1/2*A(:,end);

    Ad = A(s | d,:);
    dd = double(d(s | d));
    ss = wts(s | d).*double(s(s | d));

    % use cvx to do the constrained optimization
    if quiet 
        cvx_begin quiet
          variable delta(1) 
          variable x(Nt/2+1)
          minimize( delta )
          subject to
            -delta*dd <= Ad*x - dd <= delta*dd + delta*d2/d1*ss            
        cvx_end
    else
        cvx_begin 
          variable delta(1) 
          variable x(Nt/2+1)
          minimize( delta )
          subject to
            -delta*dd <= Ad*x - dd <= delta*dd + delta*d2/d1*ss
        cvx_end
    end
    % stack two halves together to get full linear-phase filter
    x = [x;x(end-1:-1:1)]';

%     b = x./max(abs(fft(x,osfact*Nt))); % normalized beta coefficients

    b = x./max(abs(fft(x,osfact*Nt))); % normalized beta coefficients
    
end

