%% Script for computing relations between R, P and I
% I -> the model is ignorant
% R -> the network is right
% P -> the examples is in very different perspective
% We compute the relation between the three in the cross validation set.
%
%

%% Prepare enviroment
%
clear
close all
rng(3) 
path(path,'MatFunc/Vision')
path(path,'MatFunc/Accuracy')
path(path,'MatFunc/Misc')

%% Read model data

W = load('Model/Wmom.mat').W;
U = load('Model/Umom.mat').U;
Th = load('Model/Thmom.mat').Th;
Thcn = load('Model/Thcnmom.mat').Thcn;
dict = load('Model/dictmom.mat').dict;

%% Read cross validation examples

L = size(W,2);
K = 3;
signa = size(W,1)-9;
Figures = {'Two','Three','Seven','Square','Semicircle','Star','LetterG','LetterH','LetterK'};
FLDR = 'Images/MomCV4';
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

%% Compute misclassification and certainity

error = zeros(1,Lex);
ignor = zeros(1,Lex);
for i=1:Lex
    [pred,cert] = predictcon3(W,U,s2(:,i),Th,Thcn,dict);
    if pred ~= concpt(i) 
      error(i) = 1;
    end
    ignor(i) = (cert == 0);
end

prec = 1 - sum(error)/Lex;
fprintf("The precision is: %f\n",prec);

%% Compute relations

% P(R|~I)
cpos = sum(~error & ~ignor)/sum(~ignor);
fprintf('A %f of examples that do not generate ignorance are correctly classified\n',cpos);

% perspectives annotated manually
pers = [0 0 1 1  1 1 0 0  0 1 0 1  1 1 1 0  1 1 1 0  0 1 0 1  0 1 0 1  ...
       1 0 1 1  0 1 0 1];

% P(P|~R)
cpos = sum(pers & error)/sum(error);  
fprintf('A %f of wrongly classified examples are different in perspective\n',cpos);

% P(R|~P)
cpos = sum(~pers & ~error)/sum(~pers);  
fprintf('A %f of equal in perspective examples are correctly classified\n',cpos);

% P(~I|~P) 
cpos = sum(~pers & ~ignor)/sum(~pers);  
fprintf('A %f of equal in perspective examples do not generate ignorance\n',cpos);

%% Histogram of certainty
figure;
set(gcf,'color','w');
histogram(certain,150);
xline(0.3,'--');
xlabel('Certainty');
ylabel('Counts');
title('Histogram of certainty of the model for cross validation set');