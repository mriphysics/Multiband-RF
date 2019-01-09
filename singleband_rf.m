function [rf,tb,rf_me_exc] = singleband_rf(Nt,tb,flip,mode,type,phase,d1,d2,quiet)

% if nargout <2
%     rf_me_exc = []
% end;
if and(and(strcmp(mode,'cvx'),strcmp(phase,'linear')),mod(Nt,2));
    error('Linear Design needs even number of time-points');
end
if and(strcmp(type,'exc'),flip>pi/2)
%     warning('excitation flip angle over 90 deg. Setting pulse type to "refocusing"');
%     type = 'ref';

    warning('Excitation flip angle over 90 deg. Setting flip to 90');
    flip = pi/2;
end
% Change d1 & d2 depending on pulse type ('exc','ref','me')
switch type
    case 'exc' %<-- SLR excitation pulse
        d1 = sqrt(d1/2);
        d2 = d2/sqrt(2);
    case 'ref' %<-- SLR spin-echo refocusing pulse
        d1 = d1/4;
        d2 = sqrt(d2);
    case 'me' %<-- SLR matched-exictation spin-echo refocusing pulse       
        d1 = d1/4;
        d2 = (d2/sqrt(2))^0.25;
    otherwise
        error('Pulse type not recognized');
end


% ------------------------------- %
% Change d1,d2 depending on phase arrangement ('linear','minimum','maximum','quadratic')
% Also set tbmin, if not linear.
switch phase
    case 'linear'
        fprintf('Linear phase d1=%.4f d2=%.4f\n',d1,d2);
    case {'minimum','maximum','quadratic'}
        fprintf('Non-linear phase d1=%.4f d2=%.4f\n',d1,d2);
    otherwise 
        
        error('Phase type not recognized \n Use linear, minimum, maximum or quadratic');
end
% ------------------------------- %
switch mode
    case 'cvx'
        b = dz_cvx(Nt,tb,d1,d2,quiet,phase);
    case 'ls'
        b = dz_ls(Nt ,tb,d1,d2,phase);
    otherwise
        error('Optimization mode not recognized. Choose cvx or ls');
end
    
switch phase
    case 'linear'
        
    case 'minimum'

        dinfmin = 1/2*dinf(2*d1,d2^2/2); 
        dinflin = dinf(d1,d2);      

        tbmin = tb/dinflin*dinfmin; % scale TBW product so as to get the same transition 
                                    % width as linear phase pulse with same ripples, 
                                    % after scaling back to desired slice thickness. This 
                                    % makes comparison to other MB excitations more 
                                    % meaningful, since all will have same slice characteristics.
                                    
        tb = tbmin;
        
    case 'maximum'
        
        dinfmin = 1/2*dinf(2*d1,d2^2/2); 
        dinflin = dinf(d1,d2);      

        tbmin = tb/dinflin*dinfmin; % scale TBW product so as to get the same transition 
                                    % width as linear phase pulse with same ripples, 
                                    % after scaling back to desired slice thickness. This 
                                    % makes comparison to other MB excitations more 
                                    % meaningful, since all will have same slice characteristics.
                                    
        tb = tbmin;
        
        % Time-reverse from minimum-phase to maximum-phase
        % (or equivalently, conjugate in frequency domain).
        b = b(end:-1:1);
    case 'quadratic'
        
        dinfmin = 1/2*dinf(2*d1,d2^2/2); 
        dinflin = dinf(d1,d2);      

        tbmin = tb/dinflin*dinfmin; % scale TBW product so as to get the same transition 
                                    % width as linear phase pulse with same ripples, 
                                    % after scaling back to desired slice thickness. This 
                                    % makes comparison to other MB excitations more 
                                    % meaningful, since all will have same slice characteristics.
                                    
        tb = tbmin;
        
        % Design Quadratic-phase pulse as described in Shinnar MRM 1994.
        % Evaluate roots from minimum-phase filter.
        r =roots(b);                
        % For each root with a negative imaginary component, replace it by
        % its conjugate reciprocal (i.e. "Flipping").
        r(imag(r)>=0) = 1./conj(r(imag(r) >= 0));
        % Convert from root-format to polynomial format
        b = poly(leja(r));
        % Normalise filter response.
        b = b/max(abs(freqz(b)));
end

% Apply flip angle
b = b*sin(flip/2); % scale to target flip angle    
% b = b*sin(flip/2 + atan(d1)/2); % scale to target flip angle    
% b = b*sin(flip/2 + atan(2*d1)); % scale to target flip angle    
% b = b*sin(flip/2 + atan(2*d1)); % scale to target flip angle    

[~,a] = b2amp(b,64);

% rf = -imag(islr(a,b));
rf = 1i*(islr(a,b)); %<-- need both real and imag for quadratic pulses.

if and(strcmp(type,'me'),nargout>2)
%     [~,blin] = rf2ab_sas(rf,(-N/2:N/2)-1/2,0);
    fprintf('Designing matched excitation pulse\n');
    [~,b90] = rf2ab_sas(rf,(-Nt/2:1/2:Nt/2-1/2)',0); % get beta of 180 pulse
    b90d = (b90.^2)/sqrt(2); % target 90-deg beta profile
    b90d= -conj(b90d); %negative conjugate to cancel out phase of the 180 pulse
    bx=fftshift( fft(ifftshift(b90d))/length(b90d));
    [~,ax]=b2amp(bx);
    rf_me_exc = -1i*conj(islr(ax,bx));
elseif and(~strcmp(type,'me'),nargout>2) %<--
    warning('Matched excitation pulse requested for non matched-design. Returning halve-amplitude RF')
    rf_me_exc = rf/2;
else
    rf_me_exc = [];
end

if 0
   fh = figure;
   plot(abs(fftshift(fft(b))));
   hold on;
   bx = (fftshift(fft(b)));ax = (fftshift(fft(a)));
   plot(abs(2*conj(ax).*bx));   
   legend('fft(b)','2ab')
   keyboard;
end
end

