function [ idz,pbidc ] = idmxy(z,slthick,nb,bandsep)
% Finds the ideal slice profile for a multiband pulse given a space vector
% z. getfs is too complicated, and there's no need to use frequency
% notation.

% I've made it work correctly for odd-length z, with having rndfnc = fix,
% and tps(i):tps(i+1)+1 in line 24. 
ns = length(z);
dz = z(2)-z(1);
zi = z/dz;
spos = (1:nb)-(nb+1)/2; % this is the slice shifts in units of one gap

% Passband centers
pbcs = spos*bandsep*slthick/dz*(ns-1)/ns;

% Transition points:
tps = [];
rndfnc = @(x) fix(x); %<-- define a rounding function. This definition works perfectly for odd spatial points.

pbl = floor(slthick/dz);
for i=1:length(pbcs)
	tps = [tps rndfnc(pbcs(i)-(slthick/(2*dz))+ns/2) rndfnc(pbcs(i)+slthick/(2*dz)+ns/2)];
end
idz = zeros(size(z));
j = 1;

% 08/01/2016 sas added extra code and checks to deal with the problem of
% having different passband-indices lengths. In the following for-loop a
% check is run for such case using the all_pbl_same boolean.
all_pbl_same = 1; %<-- assume all passband have equal length
for i = 1:2:2*nb;
    passband_indices = tps(i):tps(i+1)+1;
%     idz(passband_indices(1:pbl+1)) = 1;
    idz(passband_indices) = 1;
    pbidc{j} = passband_indices; %<-- first pbidc is a 1xnb cell, then it gets reshaped to a nb x pbl_max matrix.
    
    if i == 1
        first_pbl = length(passband_indices);
    elseif length(passband_indices) ~= first_pbl
        all_pbl_same = 0;
        disp('Unequal pbl found');
    end
    
    j = j+1;
end

% Find the maximum cell length
cellLength = @(x) length(x);
pbl_max = max(cellfun(cellLength,pbidc));
% When unequal passband lengths appear, repeat the final element of the
% shorter ones until they are all equal to the maximum length passband.
if all_pbl_same == 0
    for i = 1:nb
        if length(pbidc{i}) < pbl_max
            pbidc{i} = [pbidc{i}(:)' pbidc{i}(end)*ones(1,pbl_max - length(pbidc{i}))];
        end
    end
end

pbidc = reshape(cell2mat(pbidc),[pbl_max nb])';

end

