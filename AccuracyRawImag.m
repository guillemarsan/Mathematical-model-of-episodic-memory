%% Script for testing association: raw images
%

%% Prepare enviroment
%
clear
close all
rng(3)
path(path,'Images')
path(path,'MatFunc\Simulate')
path(path,'MatFunc\Patterns')
path(path,'MatFunc\Accuracy')
path(path,'MatFunc\Vision')

n = 20*20;           % neuron dimension
A = 6;           % number of neurons in concept layer
M = 200;           % number of neurons in selective layer
L = 24;          % number of stimuli

p = zeros(20,20,L);
s = zeros(n,L);
Figures = ["Triangle", "Square", "Circle", "One", "Two", "Three"];
idx = 1;
% Plot patterns
clf
for i=1:6
    for j=1:4
        Img = imread(strcat(Figures(i),int2str(j),'.jpeg'));
        p(:,:,idx) = imbox(Img,0); % center in the image
        subplot(L/4,4,idx)
        showPattern(p(:,:,idx));
        p(:,:,idx) = p(:,:,idx)/norm(p(:,:,idx)); % normalize
        aux = p(:,:,idx)';
        s(:,idx) = aux(:); % linearize
        idx = idx + 1;
    end
end
sgtitle("Images used");

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

%% Do simulations with Option D. Selective layer
%
h = 0.005;        % time step
d = 150;            % inhibitory coupling

% Sensory layer
W = SimulateNeurons4(Tmax, h, W0, s, f, alpha, b2, Th, d);

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
y = [y; y; y; y];

%U0 = 2*rand(M,A) - 1;  % random neurons
%[~,id] = sort(sum(y'*U0 > Th)); % sort neurons for convenience
%U0 = U0(:,id);

%d = 0; %no inhibition

% Concept layer
%U = SimulateNeurons3(Tmax, h, U0, y, g, alpha, bcn2, Thcn, d);

% Manually set the values
U = zeros(M*4,A);
Thcn = zeros(A,1);
for i=0:A-1
    aux = sum(y(:,i*K+1:i*K+K),2);
    U(:,i+1) = aux/norm(aux);
    v = U(:,i+1)'*y(:,i*K+1:i*K+K);
    Thcn(i+1) = min(v(v~=0));
end

%% Test accuracy

p = zeros(20,20,5);
s2 = zeros(n,5);
Figures2 = ["Triangle5", "Square5", "Circle5", "One5", "Two5", "Three5"];
figure;
for i=1:6
    Img = imread(strcat(Figures2(i),'.jpeg'));
    p(:,:,i) = imbox(Img,0);
    subplot(3,2,i)
    showPattern(p(:,:,i));
    p(:,:,i) = p(:,:,i)/norm(p(:,:,i)); % normalize
    aux = p(:,:,i)';
    s2(:,i) = aux(:); % linearize
end
V2 = W'*s2;
y2 = max(0,V2 - Th);

Vc2 = U'*[y2;y2;y2;y2]
F3 = Vc2 >= Thcn

%% Plot concept raster

V = U'*y;
F = V >= Thcn;
figure
spy(F');
title("Rasterplot concept layer neurons and stimuli they respond to");
xlabel("Neurons");
ylabel("Stimuli");