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
FLDR = 'Images/MomTrain4';
PlotFLG = true; 

figure('color','w','position',[100 100 900 900])
mom = ImportImagesEvalMoments(FLDR, Figures, signa, PlotFLG);


%% Sensory stimuli (and angles between them)

s = mom;
[n,L] = size(s);
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


%% Read crossvalidation examples

Figures = {'Two','Three','Seven','Square','Semicircle','Star','LetterG','LetterH','LetterK'};
FLDR = 'Images/MomCV4';
PlotFLG = true; 
K = 3;

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

%% Do simulations with Option D.

M = 300;          % number of neurons in the selective layer
Tmax = 300;       % max integration time
Th = 0.8;           % selective threshold
h = 0.0025;        % time step (better to decrease)
d = 150;          % inhibitory coupling
[n,L] = size(s);  % dimension and number of stimuli
f = @(t) mod(fix(t),L)+1;   % function defining the stimulus sequence
alpha = 30;  

psl = 0.975;       % selective probability
b2 = 0.7;

% Set and train the sensory layer 
W0 = 2*rand(n,M) - 1;  % random neurons
[~,id] = sort(sum(s'*W0 > Th)); % sort neurons for convenience
W0 = W0(:,id);

A = 150;       % number of neurons in the selective layer
K = 3;          % integration
Thcn = 0.125;     % conceptual threshold

pcn = 0.975; % selective probability
bcn2 = 0.25;

U0 = 2*rand(M,A) - 1;  % random neurons

pts = 5;
Thcnarr = linspace(0,0.5,pts);
bcn2arr = linspace(0,0.5,pts);
prec = zeros(pts,pts);
conscore = zeros(pts,pts);
for i = 1:pts
    for j = 1:pts
        % Hebbian learning
        [W,U] = SimulateNeurons4Layers(Tmax, h, W0, U0, s, f, K, alpha, b2, ...
                                        bcn2arr(i), Th, Thcnarr(j), d);
        y = max(0,W'*s - Th); % compute reaction to s
        V = U'*y;
        F = V >= Thcnarr(j);
        dict = conceptmap(F,K);

        error = 0;
        for k=1:Lex
            pred = predictcon4(W,U,s2(:,k),Th,Thcnarr(j),dict);
            if pred ~= concpt(k)
              error = error + 1;
            end
        end
        
        conscore(i,j) = conceptscore(dict,K);
        prec(i,j) = 1 - error/Lex;
    end
end

%% Plot results in table
fprintf('The precision is\n');
disp(prec);
fprintf('The conceptual score is\n')
disp(conscore);