%% Script for testing complex pattern recognition: digits
%

%% Prepare enviroment
%
clear
close all
path(path,'MatFunc')
path(path,'MINST')

readDigits = 40;

n = 20*20;           % neuron dimension
A = 10;           % number of neurons in concept layer
M = 200;           % number of neurons in selective layer
L = readDigits;          % number of stimuli

[p,~] = readOrdMNIST(readDigits);

% Plot patterns
clf
for i=1:readDigits
   subplot(4,10,i)
   showPattern(p(:,:,i));
end
sgtitle("Patterns used");

for j=1:readDigits
    p(:,:,j) = p(:,:,j)/norm(p(:,:,j)); % normalize
    aux = p(:,:,j)';
    s(:,j) = aux(:); % linearize
end

psl = 0.95;       % selective probability
Th = sqrt(3)*0.5; % selective threshold
pcn = 0.95;       % concept probability
Thcn = 0.5;       % concept threshold

alpha = 20;  
Tmax = 400;       % max integration time

f = @(t) mod(round(t),L)+1;   % function defining the stimulus sequence

K = 4; % associated stimuli. Must be a divisor of readDigits

% function defining the consectutive signals sequence
g = @(t) mod(round(t),L)+1-mod(round(t),K):mod(round(t),L)+1; 

delta = sqrt(1 - (2*norminv(psl) / sqrt(5*n)));
b2 = (Th/delta)^2;  % beta^2

bcn2 = (Thcn*sqrt(L)*delta*K*gamma(K + 0.5) / ...
    (Th*(1-pcn)*(1-delta)*factorial(K-1)*sqrt(M)))^2; % ???

W0 = 2*rand(n,M) - 1;  % random neurons
[~,id] = sort(sum(s'*W0 > Th)); % sort neurons for convenience
W0 = W0(:,id);

%% Do simulations with Option C. Selective layer
%
h = 0.005;        % time step (better to decrease)
d = 4;            % inhibitory coupling

% Sensory layer
W = SimulateNeurons3(Tmax, h, W0, s, f, alpha, b2, Th, d);

figure;
V = W'*s;
F = V > Th;
R = orderRasterPlot(F');
spy(R);
title("Rasterplot selective layer neurons and stimuli they respond to");
xlabel("Neurons");
ylabel("Stimuli");

%% Do simulations. Concept layer

y = max(0,W'*s - Th); % compute reaction to s
% y = y.*(y > epsilon); % avoid round to zero problems

% Manually set the values
U = zeros(M,A);
Thcn = zeros(A,1);
for i=0:A-1
    aux = sum(y(:,i*K+1:i*K+K),2);
    U(:,i+1) = aux/norm(aux);
    v = U(:,i+1)'*y(:,i*K+1:i*K+K);
    Thcn(i+1) = min(v);
end

%U0 = 2*rand(M,A) - 1;  % random neurons
%[~,id] = sort(sum(y'*U0 > Th)); % sort neurons for convenience
%U0 = U0(:,id);

%d = 0; %no inhibition

% Concept layer
%U = SimulateNeurons3(Tmax, h, U0, y, g, alpha, bcn2, Thcn, d);

% Plot raster
figure;
V = U'*y;
F = V >= Thcn;
%% R = orderRasterPlot(F');
figure
spy(F');
title("Rasterplot concept layer neurons and stimuli they respond to");
xlabel("Neurons");
ylabel("Stimuli");

%% Do simulations with Option D. 
%
% d = 150;
% 
% % Option D
% 
% W = SimulateNeurons4(Tmax, h, W0, s, f, alpha, b2, Th, d);
% 
% acD = accuracy(W,s,Th);


