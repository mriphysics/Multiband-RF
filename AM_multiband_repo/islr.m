function [ B1_slr] = islr(Ap,Bp)

% 18/04/2016 SAS New inverse SLR function named islr.m. This version
% follows the notation from Pauly's 1991 paper perfectly, and does not
% follow the Le Roux convention of taking negative conjugate of the final RF
% waveform. 

% Ap: Ntx1 vector containing alpha-polynomial coefficients.
% Bp: Ntx1 vector containing beta-polynomial coefficients.

% Output:
% B1_slr: Ntx1 vector containing SLR pulse.

% Please use under MIT license (Copyright (c) 2015 mriphysics)
% Samy Abo Seada, Anthony Price, Jo Hajnal and Shaihan Malik. January 2017
% King's College London

Nt=length(Ap);
    for j=Nt:-1:1
        A0=Ap(1);
        B0=Bp(1);
        phi(j)=2*atan2(abs(B0),abs(A0));
        theta(j)=angle(-1i*B0/A0); 
        B1_slr(j)=  phi(j)*exp(1i*theta(j));

        Cj=cos(0.5*abs(B1_slr(j)));
        Sj=1i*exp(1i*angle(B1_slr(j))) *sin(0.5*abs(B1_slr(j)));

        Ap_1=Cj.*Ap+conj(Sj).*Bp;
        Bp_1=-Sj.*Ap+Cj.*Bp;

        Ap=Ap_1(1:end-1);
        Bp=Bp_1(2:end);

    end
    
end