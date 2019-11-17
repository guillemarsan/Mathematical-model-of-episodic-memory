%% Script for studying the cluster sizes with Option A to respect to M
% 
% Sept. 30, 2019. ver. 0.1
% Given the following updating rule:
% 
% Option A:
% $$y_j(t+1) = H\left[v_j(t) - \theta - d \sum_{j\ne k} y_k(t) \right]$$
% 
% We have seen in TestInhibitoryCouplings.m that the results are far from 
% satisfactory.
% *rng(2)*
% The number of lost stimuli 
% reduces from 143 (d = 0) to 117 (d = 1). 
%
% We study whether the mean size of the clusters (group of neurons that
% respond to the same stimulus) grows lineary with respect to the number of
% neurons M. Here the sample is M = [150, 200, 250, 300, 350, 400]
% We observe that with d = 0 the relationship between M and the mean size
% of the clusters is perfectly linear: scl = [0.75, 1, 1.25, 1.5, 1.75, 2] 
% (50 in M <-> 0.25 in size)
% However, with d = 0 the slope is slightly less than 1 although not 
% significatively: scl = [0.75, 1,1.27, 1.515, 2.005]

%% Set the problem parameters
%
clear
rng(2) % for reproducibility of the results

psl = 0.95;       % selective probability
n = 30;           % neuron dimension
M = [150, 200, 250, 300, 350, 400]; % number of neurons
L = 200;          % number of stimuli
Th = sqrt(3)*0.5; % threshold
alpha = 20;  
Tmax = 200;       % max integration time
h = 0.01;         % time step
d = 1;            % inhibitory coupling

% Inhibitory function
inh = @(y) d*(sum(y) - y);

delta = sqrt(1 - (2*norminv(psl) / sqrt(5*n)));
b2 = (Th/delta)^2;  % beta^2

W0 = 2*rand(n,400) - 1;  % random neurons
s = sqrt(3/n)*(2*rand(n,L) - 1); % random stimuli
[~,id] = sort(sum(s'*W0 > Th)); % sort neurons for convenience
W0 = W0(:,id);

f = @(x) mod(round(x),L)+1;   % function defining the stimulus sequence

%% Do simulation, compute mean size cluster and show results

scl = zeros(6,1);
for i = 1:6
    W = SimulateNeurons(Tmax, h, W0(:,1:M(i)), s, f, alpha, b2, Th, inh);
    Rend = s'*W > Th;  % excitation matrix at t = Tmax
    scl(i) = mean(sum(Rend,2)); % mean size cluster 
end

% Drawing
clf
plot(M,scl,'LineWidth',2)
xlabel('neurons (M)','FontSize',14); ylabel('mean cluster size','FontSize',14)
title(['d = ' num2str(d)],'FontSize',12)

