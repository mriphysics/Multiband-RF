function [ b ] = dz_ls(Nt ,tb,d1,d2,phase);

switch phase
    case 'linear'
        
        w = dinf(d1,d2)/tb; % transition width        
        f = [0 (1-w)*(tb/2) (1+w)*(tb/2) Nt/2];%*di/dilp;
        acls=[1 1 0 0];
        wls = [1 d1/d2];

%         b =firls(2*Nt-2,f/(Nt/2),acls,wls); %Use for implementation
        b =firls(Nt,f/(Nt/2),acls,wls); %Use for implementation
        
    otherwise
        N = 2*(Nt-1);
        
        dinfmin = 1/2*dinf(2*d1,d2^2/2); 
        dinflin = dinf(d1,d2);      

        tbmin = tb/dinflin*dinfmin; % scale TBW product so as to get the same transition 
                                    % width as linear phase pulse with same ripples, 
                                    % after scaling back to desired slice thickness. This 
                                    % makes comparison to other MB excitations more 
                                    % meaningful, since all will have same slice characteristics.
        w = dinfmin/tbmin; % transition width
        
        f = [0 (1-w)*(tbmin/2) (1+w)*(tbmin/2) Nt/2];
        
        acls=[1 1 0 0];
        wls = [1 d2/d1];

%         b =firls(2*Nt-2,f/(Nt/2),acls,wls); %Use for implementation
        x =firls(N,f/(Nt/2),acls,wls); %Use for implementation
        b = real(fmp2(x));

end

end

