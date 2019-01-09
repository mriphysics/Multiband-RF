function [a,b] = rf2ab_sas(B1,gamgdt,r,Anime)
% A Bloch simulator using Cayley-Klein paramteres and the hard pulse
% approximation. Code is adapted from Grissom's blochsim.m in the
% rootflipping package, but altered to return a,b in each stage BEFORE the
% phase accrual, because this preserves the ability to show a build-up of
% Mxy over time. 

% Because gyro and dt are not included in the evaluation of C and S, for
% physical simulations scale the RF input pulse by gyro*dt

%   sas 18/05/2015

%  sas - 29/05/2015 Modify such that it runs without specifying Gradient.
%  This is such that it can run with the rootflipping code.

% Please use under MIT license (Copyright (c) 2015 mriphysics)
% Samy Abo Seada, Anthony Price, Jo Hajnal and Shaihan Malik. January 2017
% King's College London

Nt=length(B1);
if nargin == 3, %gamgdt not specified
  Anime=r;
  r = gamgdt;
  gamgdt = ones(Nt,1)*2*pi/Nt;
end;

% [Nd,Nz] = size(r); % Nz: nr of space/frequency points. Nd: nr of dimensions
[Nz,Nd] = size(r); % Nz: nr of space/frequency points. Nd: nr of dimensions

if Anime==1
    a_anime=zeros(Nz,Nt);
    a_anime(:,1)=ones(Nz,1);
    b_anime=zeros(Nz,Nt);
end

% 27/8/15 Found a bug in this code. a and b were scalars in the first
% iteration but were then turned into size Nzx1.

a=ones(Nz,1);

b=zeros(Nz,1);


for j = 1:Nt
  C = cos(abs(B1(j))/2);
  S = 1i*exp(1i*angle(B1(j)))*sin(abs(B1(j))/2);
  at = a*C - b*conj(S);
  bt = a*S + b*C;
  %Note that in Anime mode, a(1,:) contain intial conditions, and has one
  %row more than the number of iterations to account fot that.
  if Anime ==1
%       First apply phase accrual
      z= exp(-1i*r*gamgdt(j,:)');
      b_tmp=bt.*z;
%       Then correct for half angles
      z = exp(1i/2*r*sum(gamgdt(1:j,:),1)');
      a_anime(:,j)=at.*z;
      b_anime(:,j)=b_tmp.*z;
  end
  
  a = at;b = bt;
  z = exp(-1i*r*gamgdt(j,:)');
  b = b.*z;
end

if Anime ==1
%     No need for correction.
    a=a_anime;
    b=b_anime;
else
%     Same as Grissoms code.
    z = exp(1i/2*r*sum(gamgdt,1)');
    a = a.*z;
    b = b.*z; 
end