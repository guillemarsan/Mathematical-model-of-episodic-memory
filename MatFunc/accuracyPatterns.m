function [acC,acD] = accuracyPatterns(k,s)
% Given $k$, we generate a set $(6k-2)$ of binary patterns on a grid of $k
% \times k$ composed of vertical, horizontal and diagonal lines. We then
% make the model learn (option C and option D) with this stimuli and test 
% its capacity of prediction. We then repeat the experiment adding two more complex
% patterns (crosses)

psl = 0.95;       % selective probability
n = k*k;           % neuron dimension
M = 100;          % number of neurons
L = 6*k-2;          % number of stimuli
Th = sqrt(3)*0.5; % threshold
alpha = 20;  
Tmax = 400;       % max integration time

f = @(x) mod(round(x),L)+1;   % function defining the stimulus sequence

delta = sqrt(1 - (2*norminv(psl) / sqrt(5*n)));
b2 = (Th/delta)^2;  % beta^2

W0 = 2*rand(n,M) - 1;  % random neurons
[~,id] = sort(sum(s'*W0 > Th)); % sort neurons for convenience
W0 = W0(:,id);

%% Do simulations with Option C. Simple patterns
%
h = 0.005;        % time step (better to decrease)
d = 4;            % inhibitory coupling

W = SimulateNeurons3(Tmax, h, W0, s, f, alpha, b2, Th, d);

acC = accuracy(W,s,Th);

%% Do simulations with Option D. Simple patterns
%
d = 150;

% Option D

W = SimulateNeurons4(Tmax, h, W0, s, f, alpha, b2, Th, d);

acD = accuracy(W,s,Th);

end

