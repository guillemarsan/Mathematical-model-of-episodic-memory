%% Script for first steps towards conceptual layer: digits
%
% Jan. 30, 2020 ver 0.1
%   
% We use as stimuli the digits of MINST and see the results of conceptual
% learning with manually set weights for Option D:
%
% Even though the conceptual layer learns very well the different stimuli,
% it fails at classifying new examples. This leads to two conclusions:
%
% 1) Make digits more similar between them and a bigger dataset
% 2) Study lowering Tchn to avoid overfit to training examples
%
%% Prepare enviroment
%
clear
close all
path(path,'MatFunc\Simulate')
path(path,'MatFunc\Patterns')
path(path,'MatFunc\Accuracy')
path(path,'MINST')

readDigits = 20;

n = 20*20;           % neuron dimension
A = 10;           % number of neurons in concept layer
M = 400;           % number of neurons in selective layer
L = readDigits;          % number of stimuli

[p,~] = readOrdMNIST(readDigits);

% Plot patterns
clf
for i=1:readDigits
   subplot(4,10,i)
   showPattern(p(:,:,i));
end
sgtitle("Patterns used");

s = zeros(n,L);
for j=1:L
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

K = 2; % associated stimuli. Must be a divisor of readDigits

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
h = 0.005;        % time step 
d = 150;            % inhibitory coupling

% Sensory layer
W = SimulateNeurons4(Tmax, h, W0, s, f, alpha, b2, Th, d);

%% Plot selective layer

figure;
V = W'*s;
F = V > Th;
R = orderRasterPlot(F');
spy(R);
title("Rasterplot selective layer neurons and stimuli they respond to");
xlabel("Neurons");
ylabel("Stimuli");

%% Do simulations. Concept layer manually set

y = max(0,W'*s - Th); % compute reaction to s
y = [y; y; y; y];
% y = y.*(y > epsilon); % avoid round to zero problems

% Manually set the values
U = zeros(M*4,A);
Thcn = zeros(A,1);
for i=0:A-1
    aux = sum(y(:,i*K+1:i*K+K),2);
    U(:,i+1) = aux/norm(aux);
    v = U(:,i+1)'*y(:,i*K+1:i*K+K);
    Thcn(i+1) = min(v(v~=0));
end

%% Plot concept layer
figure
V = U'*y;
F = V >= Thcn;
spy(F');
title("Rasterplot concept layer neurons and stimuli they respond to");
xlabel("Neurons");
ylabel("Stimuli");

%% Test examples
[p, ~] = readMNIST("train-images.idx3-ubyte","train-labels.idx1-ubyte", 4, 370);
s2 = zeros(400,5);
for j=1:4
    p(:,:,j) = p(:,:,j)/norm(p(:,:,j)); % normalize
    aux = p(:,:,j)';
    s2(:,j) = aux(:); % linearize
end

figure()
for i=1:4
   subplot(4,10,i)
   showPattern(p(:,:,i));
end
sgtitle("Test patterns used")

V2 = W'*s2;
y2 = max(0,V2 - Th);
V2cn = U'*[y2;y2;y2;y2];
F2cn = V2cn >= Thcn; 
disp(V2cn) % output reaction to test examples
disp(F2cn)


