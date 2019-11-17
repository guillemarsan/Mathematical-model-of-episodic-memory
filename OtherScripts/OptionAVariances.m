%% Script for studying the variance of responses in a cluster with Option A
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
% We study if, after learning, the responses of all neurons active to
% a given stimulus (a cluster) are similar. 
% We observe that the mean of the variances (one variance per cluster) decreases
% with time of learning (d = 1): 3.8630e-04 (Tmax = 0), 3.9612e-07 (Tmax = 200),
% 3.8071e-12 (Tmax = 600)
% Also we see this mean is lower with d = 0: 3.8630e-04 (Tmax = 0),
% 4.4438e-10 (Tmax = 200), 2.4204e-20 (Tmax = 600)
%

%% Set the problem parameters
%
clear
rng(2) % for reproducibility of the results

psl = 0.95;       % selective probability
n = 30;           % neuron dimension
M = 300;          % number of neurons
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

W0 = 2*rand(n,M) - 1;  % random neurons
s = sqrt(3/n)*(2*rand(n,L) - 1); % random stimuli
[~,id] = sort(sum(s'*W0 > Th)); % sort neurons for convenience
W0 = W0(:,id);

f = @(x) mod(round(x),L)+1;   % function defining the stimulus sequence

%% Do simulation, compute variances and show results

W = SimulateNeurons(Tmax, h, W0, s, f, alpha, b2, Th, inh);

Y0 = max(0,s'*W0 - Th); % excitation matrix at t = 0
Yend = max(0,s'*W - Th);  % excitation matrix at t = Tmax

V0 = zeros(L,1); % initialize variance vector at t = 0
V = zeros(L,1); % initialize variance vector at t = Tmax
for j = 1:1
    vec = Y0(j,:); % do the variance for each stimulus
    V0(j) = var(vec(find(vec))); % take only the nonzero elements into account
    
    vec = Yend(j,:); % do the variance for each stimulus
    V(j) = var(vec(find(vec))); % take only the nonzero elements into account
end

mV0 = mean(V0) % compute the mean of the variances over the stimuli
mV = mean(V) 
