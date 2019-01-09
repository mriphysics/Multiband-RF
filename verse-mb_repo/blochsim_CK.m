%%% SJM 2-4-13: Bloch Simulator using CK parameters for multiple voxels
%%% simulataneously. Can also simulate single vox as long as arrays have
%%% singleton dimensions correctly defined
%%% Adapted from code by Will Grissom (http://www.vuiis.vanderbilt.edu/~grissowa/)
%%% 16-12-15: Updated to remove the '-conj(beta)' convention. Tested
%%% 16/12/15 and shown to be equivalent, except now notation is consistent
%%% with 1991 paper
%
%   [m,mz,mxyt,mzt,a,b] = blochsim_CK(B1,G,pos,sens,B0,varargin)
%
%   B1  = B1(t,c) = Nt x Nc array
%   G   = [Gx(t) Gy(t) Gz(t)] = Nt x 3 array
%   pos = [x(:) y(:) z(:)] = Ns x 3 array
%   sens= [S1(:) S2(:) ...Sn(:)] = Ns x Nc array
%   B0  = B0(:) = Ns x 1 vector
%
%   Nt = #timepoints, Ns = #space, Nc = #coils
%   B1/B0 are in mT, G in mT/m, pos in m
%
%   args: can include new parameters
%         'dt' followed by new time step in seconds (default 6.4e-6)
%         'M0' followed by new M0 (default [0;0;1])
%
%   M0 can have dimension 3xNs or 3x1. If 3xNs then pos, sens and B0 must
%   also have the correct dimensions

function [mxy,mz,mxyt,mzt,a,b] = blochsim_CK(B1,G,pos,sens,B0,varargin)

gam = 267522.1199722082;        % radians per sec per mT
dt = 6.4e-6;
M0=[0;0;1];

%%% check args
for ii=1:length(varargin)
    if strcmpi(varargin{ii},'dt')
        dt=varargin{ii+1};
    end
    if strcmpi(varargin{ii},'M0')
        M0=varargin{ii+1};
    end
end

Ns = size(pos,1);
Nt = size(G,1);

statea = ones(Ns,1);
stateb = zeros(Ns,1);

% sum up RF over coils
bxy = sens * B1.';


% sum up gradient over channels
bz = pos * G';

% add off-resonance
bz = bz + repmat(B0(:),1,Nt);


tmp = zeros(2*Ns,1);

if nargout>2
    returnallstate = 1;
    a = zeros(Ns,Nt);
    b = zeros(Ns,Nt);
else
    returnallstate = 0;
end


%%% Compute these out of loop
Phi = dt*gam*(abs(bxy).^2+bz.^2).^0.5;
Normfact = dt*gam*(Phi.^(-1));Normfact(~isfinite(Normfact)) = 0;
  
%%% now loop over time
for tt = 1:Nt
  
    phi = -Phi(:,tt); % sign reverse to define clockwise rotation
    normfact = Normfact(:,tt);
    
    nxy = normfact.*bxy(:,tt);nxy(~isfinite(nxy)) = 0;
    nz = normfact.*bz(:,tt);nz(~isfinite(nz)) = 0;
    
    cp = cos(phi/2);sp = sin(phi/2);
    alpha = cp-1i*nz.*sp;
    beta = -1i*(nxy).*sp;
    
    tmpa = alpha.*statea - conj(beta).*stateb;
    tmpb = beta.*statea + conj(alpha).*stateb;
    
    statea = tmpa;stateb = tmpb;
    if any(~isfinite(statea)); keyboard; end;
    if any(~isfinite(stateb)); keyboard; end;
    
    
    if returnallstate
        a(:,tt) = statea;
        b(:,tt) = stateb;
    end
end

% return final alpha, beta if not returning the whole 
% state progression
if ~returnallstate
  a = statea;
  b = stateb;
end

% calculate final magnetization state (M0 can be 3x1 or 3xNs)
mxy0 = M0(1,:)+1i*M0(2,:);
mz0 = M0(3,:);
mxy0=mxy0(:);
mz0 = mz0(:);

mxy = 2*mz0.*conj(statea).*stateb + mxy0.*conj(statea).^2 - conj(mxy0).*stateb.^2;
mz = mz0.*(statea.*conj(statea) - stateb.*conj(stateb));
mz = mz + 2*real(mxy0.*conj(statea).*-conj(stateb));

%%% 23-8-13: If returning all, then convert to Mxy(t) Mz(t). Assume that M0
%%% is 3x1
if returnallstate
    
    mxy0 = M0(1)+1i*M0(2);
    mz0 = M0(3);
    mxy0=mxy0(:);
    mz0 = mz0(:);
    
    mxyt = 2*mz0.*conj(a).*b;
    mxyt = mxyt + mxy0.*conj(a).^2;
    mxyt = mxyt - conj(mxy0).*b.^2;
    mzt = mz0.*(a.*conj(a) - b.*conj(b));
    mzt = mzt + 2*real(mxy0.*conj(a).*-conj(b));
end
    
    
    
end
