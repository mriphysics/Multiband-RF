function nPb= nrPb_fn(N,wp,bp)

% Code taken from Sharma et al. root-flipping code.

% Please use under MIT license (Copyright (c) 2015 mriphysics)
% Samy Abo Seada, Anthony Price, Jo Hajnal and Shaihan Malik. January 2017
% King's College London

% Function to return number of passbands in root-flipping algorithm
w = (-N/2+1:N/2)/N*2*pi;
idxPass = [];
for ii = 1:length(bp)
    idxPass = [idxPass find(w >= (bp(ii)-wp) & w <= (bp(ii)+wp))];
end
nPb=length(idxPass);