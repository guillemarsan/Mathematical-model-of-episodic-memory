%% Script for testing association: image moments
%
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

signa = 400;   % radius signature steps (moments)
Figures = {'Two','Three','Seven','Square','Semicircle','Star','LetterG','LetterH','LetterK'};
FLDR = 'Images/MomTrain3';
PlotFLG = true; 

figure('color','w','position',[100 100 900 900])
mom = ImportImagesEvalMoments(FLDR, Figures, signa, PlotFLG);


%% Sensory stimuli (and angles between them)

s = mom;
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
[n,L] = size(mom);  % dimension and number of stimuli
f = @(t) mod(round(2*t),L)+1;   % function defining the stimulus sequence
alpha = 20;  

psl = 0.975;       % selective probability
b2 = 0.7;


% Set and train the sensory layer with locality
W0 = 2*rand(n,M) - 1;  % random neurons
[~,id] = sort(sum(s'*W0 > Th)); % sort neurons for convenience
W0 = W0(:,id);

W = SimulateNeurons4(Tmax, h, W0, s, f, alpha, b2, Th, d);



%% Plot selective layer

figure('color','w','position',[100 100 1000 600])
Resp = s'*W > Th;
Resp0 = s'*W0 > Th;
PlotResultsOfSelectiveStratum(s, Resp0, Resp)

figure;
R = orderRasterPlot(Resp);
spy(R);
daspect([10 1 100]);
title("Rasterplot selective layer neurons and stimuli they respond to");
xlabel("Neurons");
ylabel("Stimuli");

%% Do simulations. Concept layer

rng(3)
A = 150;       % number of neurons in the selective layer
K = 3;          % integration
Thcn = 0.125;     % conceptual threshold
% function defining the consectutive signals sequence
g = @(t) mod(round(t),L)+1-mod(round(t),K):mod(round(t),L)+1; 
alpha = 20; 

pcn = 0.975; % selective probability
bcn2 = 0.25;

y = max(0,W'*s - Th); % compute reaction to s

% Hebbian learning
U0 = 2*rand(M,A) - 1;  % random neurons
[~,id] = sort(sum(y'*U0 > Thcn)); % sort neurons for convenience
U0 = U0(:,id);

d = 150; % no inhibition

% Concept layer
U = SimulateNeurons4(Tmax, h, U0, y, g, alpha, bcn2, Thcn, d);


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

%% Generate concept map

dict = conceptmap(F,K);

%% Read test examples

Figures = {'Two','Three','Seven','Square','Semicircle','Star','LetterG','LetterH','LetterK'};
FLDR = 'Images/Test';
PlotFLG = true; 

figure('color','w','position',[100 100 900 900])
[mom, class] = ImportImagesEvalMoments(FLDR, Figures, signa, PlotFLG);

% Change classes to concept tags
concpt = zeros(1,length(class));
for i = 0:(L/K)-1
    concpt(ismember(class,i*K+1:(i+1)*K)) = i+1;
end

s2 = mom;
[n,Lex] = size(s2);
s2 = sqrt(3/n)*(s2 - mean(s2))./std(s2);

%% Compute precision

error = 0;
for i=1:Lex
    pred = predictcon4(W,U,s2(:,i),Th,Thcn,dict);
    if pred ~= concpt(i)
      error = error + 1;
    end
end
prec = 1 - error/Lex;
fprintf("The precision is: %f\n",prec);


