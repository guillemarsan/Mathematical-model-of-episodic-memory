%% Script for testing association: raw images
%

%% Prepare enviroment
%
clear
close all
rng(3) 
path(path,'MatFunc/Vision')
path(path,'MatFunc/Accuracy')
path(path,'MatFunc/Simulate')
path(path,'MatFunc/Misc')

%% Read training data

Figures = {'One','Two','Three','Four','Five'};
FLDR = 'Images/RawTrain';
PlotFLG = true; 

figure('color','w','position',[100 100 900 900])
[imgs, class] = ImportImagesEvalRaw(FLDR, Figures, PlotFLG);


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
b2 = 0.7;


% Set and train the sensory layer with locality
W0 = 2*rand(n,M) - 1;  % random neurons
[~,id] = sort(sum(s'*W0 > Th)); % sort neurons for convenience
W0 = W0(:,id);

W = SimulateNeurons4Loc(Tmax, h, W0, s, f, alpha, b2, Th, d, loc);

figure('color','w','position',[100 100 1000 600])
PlotResultsOfSelectiveStratum(s, W0, W, Th)

%% Plot selective layer
figure;
V = W'*s;
F = V > Th;
R = orderRasterPlot(F');
spy(R);
daspect([10 1 100]);
title("Rasterplot selective layer neurons and stimuli they respond to");
xlabel("Neurons");
ylabel("Stimuli");

%% Do simulations. Concept layer

A = 600;       % number of neurons in the selective layer
K = 2;          % integration
Thcn = 0.5;     % conceptual threshold
% function defining the consectutive signals sequence
g = @(t) mod(round(t),L)+1-mod(round(t),K):mod(round(t),L)+1; 
alpha = 20; 

pcn = 0.95; % conceptual probability
bcn2 = 0.7;

y = max(0,W'*s - Th); % compute reaction to s

% Hebbian learning
U0 = 2*rand(M,A) - 1;  % random neurons
[~,id] = sort(sum(y'*U0 > Thcn)); % sort neurons for convenience
U0 = U0(:,id);

loc = A/20;
d = 150; %no inhibition

% Concept layer
U = SimulateNeurons4Loc(Tmax, h, U0, y, g, alpha, bcn2, Thcn, d, loc);

%% Plot concept layer
V = U'*y;
F = V >= Thcn;
R = orderRasterPlot(F');
figure
spy(R);
daspect([10 1 100]);
title("Rasterplot concept layer neurons and stimuli they respond to");
xlabel("Neurons");
ylabel("Stimuli");

return
%% Read test examples

Figures = {'One','Two','Three','Four','Five'};
FLDR = 'Images/RawTest';
PlotFLG = true; 

figure('color','w','position',[100 100 900 900])
[imgs, class] = ImportImagesEvalRaw(FLDR, Figures, PlotFLG);

s2 = imgs;
[n,Lex] = size(s2);
s2 = sqrt(3/n)*(s2 - mean(s2))./std(s2);


%% Compute precision

V2 = W'*s2;
y2 = max(0,V2 - Th);

Vcn2 = U'*y2;
F3 = Vcn2 >= Thcn;
figure
disp("Reaction of concept layer to test examples:")
spy(F3')
daspect([10 1 100]);

% Compute precision
error = 0;
for i=1:Lex
    if predictcon(W,U,s2(:,i),Th,Thcn) ~= i
      error = error + 1;
    end
end
prec = 1 - error/Lex;
fprintf("The precision is: %f\n",prec);