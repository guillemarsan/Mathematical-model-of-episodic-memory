%% Script to study stimulus and selective response: raw images
%
% Script to study how good the option to use boxed images is.
% We see that pictures of the same object even with high inclination or 
% rotation are very similar between each other and not very ortogonal.
% Therefore, with this method is better to use one example of each symbol
% to train. 
% Using 5 examples of each symbol : Sel = 96 Inact = 24 Rmean = 2
% Cherry picking examples: Sel = 155 Inact = 36 Rmean = 1.2
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

Figures = {'One','Two','Three','Four','Five'};
PlotFLG = true; 

figure('color','w','position',[100 100 900 900])
[imgs, class] = ImportImagesEvalRaw(Figures, PlotFLG);


%% Sensory stimuli (and angles between them)

s = imgs;
[n,~] = size(s);
s = sqrt(3/n)*(s - mean(s))./std(s);

nrmS = sqrt(sum(s.^2)); % norma s
S = s./nrmS;
CosAngle = S'*S; % cos(angle)

figure('color','w')
imagesc(CosAngle,[-1 1])
axis square
colorbar
colormap('jet')
xlabel('stimulus','FontSize',14)
ylabel('stimulus','FontSize',14)
title('cosine of the angles among pairs of stimuli','FontSize',14)

%% Do simulations with Option D. Selective layer

M = 300;          % number of neurons in the selective layer
Tmax = 400;       % max integration time
Th = 0.8;           % selective threshold
h = 0.0025;        % time step (better to decrease)
d = 150;          % inhibitory coupling
loc = M/20;       % locality of inhibition
[n,L] = size(imgs);  % dimension and number of stimuli
f = @(t) mod(round(2*t),L)+1;   % function defining the stimulus sequence
alpha = 20;  

psl = 0.975;       % selective probability
delta = sqrt(1 - (2*norminv(psl) / sqrt(5*n)));
b2 = (Th/delta)^2;  % beta^2


% Set and train the sensory layer with locality
W0 = 2*rand(n,M) - 1;  % random neurons
[~,id] = sort(sum(s'*W0 > Th)); % sort neurons for convenience
W0 = W0(:,id);

W = SimulateNeurons4Loc(Tmax, h, W0, s, f, alpha, b2, Th, d, loc);

figure('color','w','position',[100 100 1000 600])
PlotResultsOfSelectiveStratum(s, W0, W, Th)