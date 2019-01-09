function peak = flip2peak(N,r,wp,bp,p,symtype,d1)

% 30/09/2016 sas - New flip2peak function that can deal with all three symmetry
% types. This new one simply calls flip2rf and only returns the peak such
% that it can be used as a cost function for a Genetic algorithm.

% % % % % Old description: ignore.
% 07/09/2015 sas - function that takes in as input a binary vector p of
% length, the number of passbands (or half if symmetric constraint is
% imposed) and returns a peak value. 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Adapted code from Grissoms flipZeros function, which in turn was adapted
% from Miki Lustig. This function will be used in a genetic algorithm
% implementation of root-flipping

% Please use under MIT license (Copyright (c) 2015 mriphysics)
% Samy Abo Seada, Anthony Price, Jo Hajnal and Shaihan Malik. January 2017
% King's College London

[~,peak,~] = flip2rf(N,r,wp,bp,p,symtype,d1);