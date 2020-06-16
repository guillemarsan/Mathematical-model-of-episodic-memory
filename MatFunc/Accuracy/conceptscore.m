function score = conceptscore(dict)
    % Function which, given the dictionary of conceptual neurons computes a
    % 'score' 
    %
    % Inputs
    %   dict = 1xa array where each value is the concept each conceptual
    %   neuron associates. -1 if it is not conceptual
    % Outputs
    %   score = integer between 0 and n+1 as follows:
    %       0: no concepts are learned
    %       1,2,...,n: number of concepts learned
    %       n+1: all concepts are learned and there are spare neurons for
    %       learning other new concepts
    
    concepts = unique(dict);
    A = length(dict);
    minn = A/20;
    newcon = A;
    score = 0;
    for i = 2:length(concepts)
        numcon = sum(dict == concepts(i));
        if numcon > minn
            score = score + 1;
        end
        if numcon < newcon
            newcon = numcon;
        end
    end
    if score == length(concepts)-1
        if sum(dict == -1) > newcon
            score = score+1;
        end
    end
end

