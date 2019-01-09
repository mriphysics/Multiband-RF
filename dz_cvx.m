function [b] = dz_cvx(Nt,tb,d1,d2,quiet,phase)

switch phase
    case 'linear'
    osfact = 10; % oversampling factor for frequency grid

    nn = (0:Nt/2*osfact)'/(Nt/2*osfact);  % 0 to 1 - profile indices
    d = zeros(Nt/2*osfact+1,1);          % passband mask
    s = zeros(Nt/2*osfact+1,1);          % stopband mask
    wts = zeros(Nt/2*osfact+1,1);        % ripple taper weights

    w = dinf(d1,d2)/tb;
    % start out the f, m and w vectors with the DC band
    f = [0 (1-w)*(tb/2) (1+w)*(tb/2)];%*di/dilp;
    d = nn <= f(2)/(Nt/2); % target pattern
    wts = 1./abs(nn).^2; % quadratically-decaying ripple weights
    % append the last stopband
    s = s | (nn >= f(end)/(Nt/2));
    wts = wts.*s;wts(isnan(wts)) = 0;wts = wts./max(wts);

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
    
    otherwise

%     d1 = 0.01;
%     d2 = 0.01;
    

    N = 2*(Nt-1);
    osfact = 20; % oversampling factor

    % Apply parameter relations for spin-echo:
%     d1 = d1/4;
%     d2 = sqrt(d2);

    % For a lin-phase based design, double passband ripple and square
    % passband ripple

%     d1 = 2*d1;
%     d2 = d2^2/2; %<-- makes the passband ripple the same, but double the stopband ripple (Still below
        
    nn = (0:N/2*osfact)'/(N/2*osfact);  % 0 to 1 - profile indices
    d = zeros(N/2*osfact+1,1);          % passband mask
    s = zeros(N/2*osfact+1,1);          % stopband mask
    wts = zeros(N/2*osfact+1,1);        % ripple taper weights
        
    dinfmin = 1/2*dinf(2*d1,d2^2/2); 
    dinflin = dinf(d1,d2);      

    tbmin = tb/dinflin*dinfmin; % scale TBW product so as to get the same transition 
                                % width as linear phase pulse with same ripples, 
                                % after scaling back to desired slice thickness. This 
                                % makes comparison to other MB excitations more 
                                % meaningful, since all will have same slice characteristics.
    w = dinfmin/tbmin; % transition width
    % start out the f, m and w vectors with the DC band
    f = [0 (1-w)*(tbmin/2) (1+w)*(tbmin/2)];
    
    d = nn <= f(2)/(Nt/2); % target pattern
    wts = 1./abs(nn).^2; % quadratically-decaying ripple weights
    % append the last stopband
    s = s | (nn >= f(end)/(Nt/2));
    wts = wts.*s;wts(isnan(wts)) = 0;wts = wts./max(wts);

    % build system matrix for cvx design
    A = 2*cos(2*pi*(0:N/2*osfact)'*(-N/2:0)/(N*osfact));
    A(:,end) = 1/2*A(:,end);

    Ad = A(s | d,:);
    dd = double(d(s | d));
    ss = wts(s | d).*double(s(s | d));

    % use cvx to do the constrained optimization
    if quiet == 1
        cvx_begin quiet
          variable delta(1) 
          variable x(N/2+1)
          minimize( delta )
          subject to
              -delta*dd <= Ad*x - dd <= delta*dd + delta*d2^2/(2*d1)*ss           
        cvx_end
    else        
        cvx_begin
          variable delta(1) 
          variable x(N/2+1)
          minimize( delta )
          subject to
              -delta*dd <= Ad*x - dd <= delta*dd + delta*d2^2/(2*d1)*ss            
        cvx_end
    end
    x = [x;x(end-1:-1:1)]';

    blin=x;
    % factor the linear phase filter to get a min-phase filter b    
    
    b = real(fmp2(x));
%     b = b(end:-1:1);
end

