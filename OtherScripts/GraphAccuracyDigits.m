%% Script for testing complex pattern recognition: digits
%
% Nov. 17, 2019 ver 0.1
% 
% We compute the mean and standard deviation of the accuracy of OptionC and
% OptionD over a some (readDigits) MINST database patterns (20x20
% images of digits)
%
% *Option D:*
%
% $$v_{j}^{in}(t) = d \sum_{k=1}^{M} y_k(t)(y_k(t) - y_j(t)) $$ 
%
% $$y_{j}^{in}(t) = H(v_{j}^{in}(t) - \theta_{in})$$
%
% $$v_{j}(t) = w_{j}(t)s(t) - 1 * y_{j}^{in}(t)$$
%
% $$y_j(t+1) = H\left[v_j(t) - \theta \right]$$
% 
% *Option C:*
% 
% $$\dot{w}_j = \alpha y_j\left(\left[\beta^2 - d\sum_{k\ne j} y_k\right]s - v_jw_j\right)$$
%
% $$y_j = H(v_j - \theta)$$
%
% *readDigits = 20*
%           
% We get very good results for both except the mean is greater for OptionD 
% and has a considerably smaller standard deviation:
%
% Mean accuracy OptionC: 0.795
% Standard deviation OptionC: 0.1301
%
% Mean accuracy OptionD: 0.88
% Standard deviation OptionD: 0.0823

%% Compute results
%
clear
close all
path(path,'MatFunc')
path(path,'MINST')

readDigits = 20;

n = 20*20;           % neuron dimension
M = 100;          % number of neurons in selective layer
L = readDigits;          % number of stimuli

it = 10; % number of iterations

resultsC = zeros(1,it);
resultsD = zeros(1,it);
for i = 1:it
    i
    rng(i); % use different initial values
    offset = randi(30000); % offset to read in the database
    
    % read the database
    [p, ~] = readMNIST("train-images.idx3-ubyte","train-labels.idx1-ubyte", readDigits, offset);

    for j=1:readDigits
        p(:,:,j) = p(:,:,j)/norm(p(:,:,j)); % normalize
        aux = p(:,:,j)';
        s(:,j) = aux(:); % linearize
    end

    [acC,acD] = accuracyPatterns(n,M,L,s);
    resultsC(i) = acC;
    resultsD(i) = acD;
end

meanC = mean(resultsC) % mean
stdevC = std(resultsC) % standard deviation
meanD = mean(resultsD)
stdevD = std(resultsD)
