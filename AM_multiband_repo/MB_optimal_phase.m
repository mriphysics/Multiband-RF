%%% Simple script to compute optimal slice phase offsets for multi-band
%%% imaging.
%%% Generates results presented in ISMRM 2015 abstract #2398 
%%% Shaihan Malik 12-5-2015


%% Compute optimized slice phases for specified MBF

%%%%%%%% User Inputs %%%%%%%%%%%%%%
MBF = 6;                %<--- Number of slices
am_only = true;         %<--- AM only (true) or AM/FM (false)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M=512; %<- number of time points, arbitrary
spos = (1:MBF)-(MBF+1)/2; %<- relative slice positions

%%% Function for max amplitude of arbitrary set of phases
bt = @(phi)(sum(exp(2*pi*1i*linspace(0,1,M)'*spos).*repmat(exp(1i*phi),[M 1]),2));

%%% If am_only need function to map phases to slice pairs
if am_only

    %%% identify pairs of slices that are equally offset (odd/even separate)
    if mod(MBF,2)==0 
        ixmap = @(x)(([x -fliplr(x)]));
    else
        ixmap = @(x)(([x -fliplr(x(1:end-1))]));
    end
else
    % If not mapping slices, just return input
    ixmap = @(x)(x);
end
bmax = @(phi)(max(abs(bt(ixmap(phi)))));


%%% Cost function for optimization
if am_only
    % Maximum amplitude
    costfun = @(phi)(max(abs(real(bt(ixmap(phi)))))+100*max(abs(imag(bt(ixmap(phi))))));
    nsol = ceil(MBF/2); %<-- the number of phases we need to find (searching for pairs)
else
    costfun = bmax;
    nsol = MBF;
end

%%% Run minimization multiple times with diff initial guess
nstarts=10;
cf=[];
phi_sol = zeros([1 nsol]);
opt=optimset('Display','none');
for ii=1:nstarts
    
    phi0 = rand([1 nsol]);
    
    phis = fmincon(costfun,phi0,[],[],[],[],zeros(size(phi0)),2*pi*ones(size(phi0)),[],opt);
    cf(ii) = costfun(phis);
    
    % Check if this is the lowest cost solution
    if costfun(phis)<costfun(phi_sol)
        phi_sol = phis;% this is now the lowest cost, so replace
    end
   fprintf('iteration %d, cf=%1.3f current best=%1.3f\n',ii,cf(ii),costfun(phi_sol))
end
    
%%% Report the best solution
fprintf(1,['\nMBF=%d, bmax = %1.3f, PHASES (radian): ' repmat('%1.3f ',[1 MBF]) '\n'],...
    MBF,bmax(phi_sol),ixmap(phi_sol));


%% Plot modulation function
figure;
subplot(3,1,1)
if am_only
    plot(real(bt(ixmap(phi_sol))))
else
    plot(abs(bt(ixmap(phi_sol))))
end
grid on
xlim([1 M])
ylabel('Amplitude / au')
title('Amplitude modulation function')
set(gca,'xticklabel','')

subplot(3,1,2)
if am_only
    plot(imag(bt(ixmap(phi_sol))));ylim([-1 1])
    %(this is for display purposes to avoid pi phase jumps in a real waveform. 
    % Try yourself if you think this is a cheat)
else
    plot(angle(bt(ixmap(phi_sol))))
end
grid on
xlim([1 M])
ylabel('phase / rad')
title('Phase modulation function')
set(gca,'xticklabel','')

subplot(3,1,3)
hold on
plot(cos(linspace(0,2*pi,M)),sin(linspace(0,2*pi,M)))
plot(cos(ixmap(phi_sol)),sin(ixmap(phi_sol)),'or')
grid on
title('Slice phases in complex plane')

set(gcf,'position',[100 100 250 750])
