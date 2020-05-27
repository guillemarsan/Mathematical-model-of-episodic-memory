%% Script to locality of inhibition of Option D: image moments
%
%
% We compare how different locality of inhibition results in different
% number of selective neurons. We check as well the difference between 
% having a hard boundary vs window locality.
%
% rng(3)
%
% We show that neither of them perform much better than global inhibition
%
%
%% Prepare enviroment
%
clear
close all
rng(3) 
path(path,'MatFunc/Vision')
path(path,'MatFunc/Simulate')
path(path,'MatFunc/Misc')

%% Read data

signa = 400;   % radius signature steps (moments)
FLDR = 'Images/MomTrain';
Figures = {'One','Two','Three','Four','Triangle','Square'};
PlotFLG = true; 

figure('color','w','position',[100 100 900 900])
[mom, class] = ImportImagesEvalMoments(FLDR, Figures, signa, PlotFLG);

s = mom;
[n,~] = size(s);
s = sqrt(3/n)*(s - mean(s))./std(s);

%% Do simulations with haard boundary and sliding window

M = 300;          % number of neurons in the selective layer
Tmax = 400;       % max integration time
Th = 0.8;           % selective threshold
h = 0.0025;        % time step (better to decrease)
d = 150;          % inhibitory coupling
[n,L] = size(mom);  % dimension and number of stimuli
f = @(t) mod(round(2*t),L)+1;   % function defining the stimulus sequence
alpha = 20;  

psl = 0.975;       % selective probability
delta = sqrt(1 - (2*norminv(psl) / sqrt(5*n)));
b2 = (Th/delta)^2;  % beta^2

% Set and train the sensory layer with locality
W0 = 2*rand(n,M) - 1;  % random neurons
[~,id] = sort(sum(s'*W0 > Th)); % sort neurons for convenience
W0 = W0(:,id);

karr = [300, 150, 100, 50, 25, 3, 2, 1];
i = 1;
Sel = zeros(1,length(karr));
SelHard = zeros(1,length(karr));
for k = karr
    
    loc = M/k;       % locality of inhibition
    
    % Sliding window
    W = SimulateNeurons4Loc(Tmax, h, W0, s, f, alpha, b2, Th, d, loc);
    
    y = SimulateConvergence1(W,s,Th,d,0.01,1000,loc,false);
    Resp = y > 0;
    R = sum(Resp,2);
    Sel(i) = sum(R==1);
    
    % Hard boundary
    W2 = SimulateNeurons4LocHard(Tmax, h, W0, s, f, alpha, b2, Th, d, loc);
    
    y2 = SimulateConvergence1(W2,s,Th,d,0.01,1000,loc,true);
    Resp2 = y2 > 0;
    R2 = sum(Resp2,2);
    SelHard(i) = sum(R2==1);
    i = i+1;
end

%% Do simulation without locality

W3 = SimulateNeurons4(Tmax, h, W0, s, f, alpha, b2, Th, d);
y3 = SimulateConvergence2(W3,s,Th,d,0.01,1000);
Resp3 = y3 > 0;
R3 = sum(Resp3,2);
SelNoLoc = sum(R3==1);

%% Plot results
figure;
hold on
plot(karr, ones(size(karr))*SelNoLoc,'k');
plot(karr,Sel,'r');
plot(karr,SelHard,'b');
ylim([0 300])
xlim([1 300])
title('Selectivity after learning with different locally inhibition rules')
xlabel('group / sliding window size','FontSize',16)
ylabel('# selective neurons','FontSize',16)
legend('No locality','Hard boundary','Sliding window')
