function [acC,acD] = accuracyPatterns(n,M,L,s)
% Given the stimuli s we compute the accuracy for OptionC and OptionD

psl = 0.95;       % selective probability
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

