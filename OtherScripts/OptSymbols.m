%% Script to analyse the symbols and pick the most distinguishable
%
% We study which are the best symbols to work with. We set a threshold for
% the cosine of the angle the stimuli form. If the cosine is smaller we
% consider that the algorithm is not going distinguish them. We compute the
% least error a symbol can get when trying to distinguish it from other
% symbols and when trying to associate it with other examples of the same 
% symbol. We then output the symmax symbols of each class with least sum of
% errors
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
Figures = {'One','Two','Three','Four','Five','Seven','Triangle','Square','Circle','Semicircle',...
       'Star','Bar','LetterV','LetterT','LetterK','LetterH','LetterL','LetterG'};
PlotFLG = false; 

figure('color','w','position',[100 100 900 900])
[mom, class] = ImportImagesEvalMoments(FLDR, Figures, signa, PlotFLG);


%% Algorithm

s = mom;
[n,L] = size(s);
s = sqrt(3/n)*(s - mean(s))./std(s);

nrmS = sqrt(sum(s.^2)); % norm s
S = s./nrmS;

K = 30;                 % symbols per class
sym = length(Figures);  % symbols
numex = L/sym;          % examples per symbol

thresh = 0.5;

symbol = zeros(1,sym);
intersymbol = zeros(1,sym);
j = 1;
for i = 1:numex:L
    idx = i:i+numex-1;          % indeces of examples of same symbol
    cidx = setdiff(1:L,idx);    % indecces of examples of different symbol
    
    % min number of false negatives
    angles = S(:,idx)'*S(:,idx);
    symbol(j) =  min(sum(angles<thresh))/numex; 
    
    % min number of false positives
    anglesinter = S(:,idx)'*S(:,cidx);
    intersymbol(j) = min(sum(anglesinter>thresh,2))/(L-numex);
    
    j = j+1;
end

%% Plot results

figure;
set(gcf,'color','w');
hold on;
c = repelem(1:(L/K),K/numex);
gscatter(symbol,intersymbol,c);
xlim([0,1]);
ylim([0,1]);
xlabel('Best false negative error');
ylabel('Best false positive error');
title('Distribution of symbols');
text(symbol+0.01, intersymbol+0.01, Figures, 'FontSize',8);
legend('Digits','Shapes','Letters');
plot(0:0.1:1, 0.5*ones(1,11), '--k','HandleVisibility','off');
plot(0.5*ones(1,11), 0:0.1:1, '--k','HandleVisibility','off');

%% Choose the best symmax symbols of each class

symmax = 3;
j = 1;
idx = zeros(1,symmax*(L/K));
for i = 1:(K/numex):sym
    err = symbol(i:i+(K/numex)-1) + intersymbol(i:i+(K/numex)-1);
    errs = sort(err);
    errmax = errs(1:symmax);
    idx(j:j+(symmax-1)) = find(ismember(err,errmax))+(i-1);
    j = j+symmax;
end

%% Print results

tags = Figures(idx)';
fprintf('The best symbols to use are:\n');
fprintf('%s\n',tags{:});
