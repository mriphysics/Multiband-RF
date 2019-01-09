%% 30-9-13: Distort gradients to predict effect of 43us eddy current term
%%% Function also makes transmit and receive kspaces 
% modified version: if tau=0 then return unmodified grads along with
% k-spaces (this can be used to quickly get tx/rx k-spaces)
function [Gdistorted,ktx,krx] = gradient_distort_FT(Ginput,tau,delay,dt,Npad)

if tau==0
    runDistort=false;
else
    runDistort=true;
end

if numel(tau)<3
    tau = tau(1)*[1 1 1]; % allow channel dependent tau
end
if ~exist('delay','var')
    delay=tau;
end
if ~exist('dt','var')
    dt=6.4e-6;
end
    
%%% add some zeros on front and end
if ~exist('Npad','var')
    Npad=100;
end


Ngradients = size(Ginput,2);
G = cat(1,zeros([Npad Ngradients]),Ginput,zeros([Npad Ngradients]));
M = length(G);
t = (0:M-1)*dt;

ff = -0.5/dt:1/((M-1)*dt):0.5/dt;

    % Find required frequency resolution
    df_act = ff(2)-ff(1);
    %%% construct Fourier matrix
    F = dt*exp(2*pi*1i*ff'*t);
    % Also define inverse
    Fi = df_act*exp(-2*pi*1i*t'*ff);    

    H = repmat(1./tau,[length(ff) 1]).*(F*exp(-repmat(t(:),[1 3])./repmat(tau,[M 1])));
    H = H./max(abs(H(:)));
    
    %%% Apply transfer function
    Gdistorted = zeros([M Ngradients]);
    %%% check if H has all axes
    if size(H,2)==1
        H=repmat(H,Ngradients);
    end
    if size(H,2)~=Ngradients
        fprintf(1,'Error: GIRF has different number of axes to gradient\n')
        return
    end
    
    %%% Apply H
    Gdistorted = real(Fi*(H.*(F*G)));
   
    %%% Apply Delay
    for ix = 1 : Ngradients
        Gdistorted(:,ix) = interp1(t,Gdistorted(:,ix),t+delay(ix));
    end
    
%     debug:
%     girfll = load('GIRF_all_axes_20140729');
%     figure;
%     FG = abs(F*G(:,3));
%     plot(ff,FG/max(FG));
%     hold on;
%     plot(ff,abs(H(:,3)));
%     plot(girfll.ff,abs(girfll.Hw(:,3)));
%     dum = 1;   

%%% remove padding
Gdistorted([1:Npad (M-Npad+1):M],:)=[];
M = length(Gdistorted);
t = (0:M-1)*dt;

%%% Now make k-spaces
gamma_mT = 2*pi*42577.46778; % rad s^-1 mT^-1

%%% Tx - integrate to end
ktx = -gamma_mT*flipud(cumtrapz(t,flipud(Gdistorted)));
%%% Rx - simple integration as you go along
krx = gamma_mT*cumtrapz(t,Gdistorted);

end