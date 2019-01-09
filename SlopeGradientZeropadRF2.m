function [Gs,RFs] = SlopeGradientZeropadRF2(G,RF,slewrate,dt)
% 14/03/16 This version assums that G has non-zero edges, and simply
% appends a linear gradient to the first and last value subject to a
% maximum slew rate


% Takes in a column vector of Gradient and RF waveform, as well as their
% dwell-time and the maximum allowable slew rate. 

% Output A sloped gradient subject to maximum slew-rate and a zero-padded
% RF waveform such that the two outputs have equal length.

l = length(RF);
l2 = length(G);
if l ~= l2
    error('RF and G of un-equal length')
end

% Smooth gradient
 
% Find the amount of up-hill and down-hill elements
 find(max(G(1:fix(l2/2))));
 uphill_index = ceil(G(1)/(slewrate*dt))+1;
 downhill_index = ceil(G(end)/(slewrate*dt))+1;
 
 dum = 1;

 Gs = [linspace(0,G(1),uphill_index)';G(2:end-1);linspace(G(end),0,downhill_index)'];
 % Append zeros left and right of the RF waveform
 RFs = [zeros(uphill_index-1,1);RF(:);zeros(downhill_index-1,1)];
 
waveforms = [Gs./max(Gs) RFs./max(RFs)];

end

