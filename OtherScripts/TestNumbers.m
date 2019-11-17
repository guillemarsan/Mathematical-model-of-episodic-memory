%% Script for testing complex pattern recognition: digits
%
% Nov. 17, 2019, ver 0.1
%
% We read a number of patterns from the MINST database (20x20 images of
% numbers) and test the accuracy for OptionC and OptionD at learning this
% patterns
% 
% *Option D:*
%
% $$v_{j}^{in}(t) = d \sum_{k=1}^{M} y_k(t)(y_k(t) - y_j(t)) $$ 
%
% $$y_{j}^{in}(t) = H(v_{j}^{in}(t) - \theta_{in})$$
%
% $$v_{j}(t) = w_{j}(t)s(t) - 1 * y_{j}^{in}(t)$$
%
% $$y_j(t+1) = H\left[v_j(t) - \theta \right]$$
% 
% *Option C:*
% 
% $$\dot{w}_j = \alpha y_j\left(\left[\beta^2 - d\sum_{k\ne j} y_k\right]s - v_jw_j\right)$$
%
% $$y_j = H(v_j - \theta)$$
%
% Oct. 23, 2019, ver 0.2.
%
% *rng(3)*
% *readDigits = 20*
%
% The results are very positive since OptionC gets 0.9 of accuracy and
% OptionD gets 1 of accuracy.
%

%% Set the problem parameters
%
clear
close all
path(path,'MatFunc')
path(path,'MINST')
rng(3) % for reproducibility of the results

readDigits = 20;    % number of digits read
psl = 0.95;       % selective probability
n = 20*20;           % neuron dimension
M = 200;          % number of neurons
L = readDigits;          % number of stimuli
Th = sqrt(3)*0.5; % threshold
alpha = 20;  
Tmax = 400;       % max integration time

f = @(x) mod(round(x),L)+1;   % function defining the stimulus sequence

delta = sqrt(1 - (2*norminv(psl) / sqrt(5*n)));
b2 = (Th/delta)^2;  % beta^2

% read from the database
[p, ~] = readMNIST("train-images.idx3-ubyte","train-labels.idx1-ubyte", readDigits, 0);

for i=1:readDigits
    p(:,:,i) = p(:,:,i)/norm(p(:,:,i)); % normalize
    aux = p(:,:,i)';
    s(:,i) = aux(:); % linearize
end
    
W0 = 2*rand(n,M) - 1;  % random neurons
[~,id] = sort(sum(s'*W0 > Th)); % sort neurons for convenience
W0 = W0(:,id);

% Plot some patterns
clf
for i=1:8
   subplot(4,2,i)
   showPattern(p(:,:,i*floor(L/8)));
end
sgtitle("Some of the patterns used");
drawnow 

%% Do simulations with Option C.
%
h = 0.005;       % time step (better to decrease)
d = 4;            % inhibitory coupling

W = SimulateNeurons3(Tmax, h, W0, s, f, alpha, b2, Th, d);

ac = accuracy(W,s,Th)


% Plot raster
figure;
V = W'*s;
F = V > Th;
R = orderRasterPlot(F');
spy(R);
title("Rasterplot neurons and stimuli they respond to");
xlabel("Neurons");
ylabel("Stimuli");

% Plot patterns
figure
j = 1;
for i=1:5
    subplot(5,2,j); j = j + 1;
    showPattern(p(:,:,i*floor(L/5)));
    pred = predict(V,W,s(:,i*floor(L/5)),Th);
    subplot(5,2,j); j = j + 1;
    showPattern(p(:,:,pred));
end 
sgtitle("Shown stimulus vs predicted stimulus");


%% Do simulations with Option D
%
h = 0.005;
d = 150;

% Option D

W = SimulateNeurons4(Tmax, h, W0, s, f, alpha, b2, Th, d);

ac = accuracy(W,s,Th)

% Plot raster
figure;
V = W'*s;
F = V > Th;
R = orderRasterPlot(F');
spy(R);
title("Rasterplot neurons and stimuli they respond to");
xlabel("Neurons");
ylabel("Stimuli");

% Plot patterns
figure
j = 1;
for i=1:5
    subplot(5,2,j); j = j + 1;
    showPattern(p(:,:,i*floor(L/5)));
    pred = predict(V,W,s(:,i*floor(L/5)),Th);
    subplot(5,2,j); j = j + 1;
    showPattern(p(:,:,pred));
end 
sgtitle("Shown stimulus vs predicted stimulus");
