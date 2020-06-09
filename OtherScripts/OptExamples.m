%% Script to analys the training set and pick the most distinguishable
%
% We pick greedily the examples from the training set such that the cosine
% of the angle with the already picked set is less than a certain threshold
% to ensure that they are orthogonal enough.
%
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

signa = 400;   % radius signature steps (moments)
FLDR = 'Images/Experiment';
Figures = {'Two','Three','Seven','Square','Semicircle','Star',...
            'LetterK','LetterH','LetterG'};
PlotFLG = false; 

figure('color','w','position',[100 100 900 900])
[mom, class] = ImportImagesEvalMoments(FLDR, Figures, signa, PlotFLG);


%% Algorithm
% Greedy algorithm to pick most distinguishable examples
s = mom;
[n,L] = size(s);
s = sqrt(3/n)*(s - mean(s))./std(s);

nrmS = sqrt(sum(s.^2)); % norm s
S = s./nrmS;

sym = length(Figures); % number of symbols
numex = L/sym;         % number of examples per symbol

thresh = 0.55;         % threshold of acceptance

idx = [1];
for i = (numex+1):numex:L
    for j = 0:(numex-1)
        cos = S(:,(i+j))'*S(:,idx);
        if ~(max(abs(cos)) > thresh) % if the cosine of all the angles < thresh
            idx = [idx,(i+j)];       % we add it to the set
            break
        end
    end
end

%% Print results

if length(Figures) > length(idx)
    % If we cannot find an example for each symbol the training set is
    % insufficient
    fprintf('Could not find one example for each symbol');
else
    fprintf('The best training examples are\n');
    for i = 1:length(idx)
        ex = mod(idx(i),numex);
        if ex == 0 
            ex = 5;
        end
        fprintf('%s%i\n',Figures{fix(idx(i)/numex)+1},ex);
    end
end
