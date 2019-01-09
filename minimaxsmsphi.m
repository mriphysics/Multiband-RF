function [ phir,maxerror ] = minimaxsmsphi( mb,idc,Mxy,z)
% 11/08/2016 
% A function that, given a multiband factor, passband indices, slice
% profile and space vector z (on which Mxy is defined) finds the rewind
% area to minimise the maximum phase difference across all the slices

% sas

% 22/03/18 Minimize std directly.

cost = @(phir) nestfun(phir);
phir0 = 1;
phir = fminsearch(cost,phir0);

if 0
    figure;
    plot(z,abs(Mxy),'b');
    hold on;
    scatter(z(idc),abs(Mxy(idc)),'r');
    dum = 1;
end

maxerror = nestfun(phir);
function maxerr = nestfun(phir)
    maxerr = 0;
    for i = 1:mb
        psb = idc(2*i-1):idc(2*i);
        phs = unwrap(angle(Mxy(psb).* exp(-1i*phir*z(psb))));
%         err = max ( abs( phs - mean(phs) ) );
% 22/03/18 Minimize std directly.
        err = max ( std(phs) );
        if err > maxerr
            maxerr = err;
        end
    end
end

end