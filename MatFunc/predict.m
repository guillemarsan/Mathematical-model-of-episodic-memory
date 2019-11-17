function pred = predict(V,W,s,Th);
    % Function which, given a stimulus s, predicts 
    % the number of the stimulus
    %
    % Inputs
    %   V = MxL matrix representing the reponse of 
    %   each neuron to each stimulus
    %   W = nxM weights matrix
    %   s = nx1 stimulus shown
    %   Th = theta parameter
    % Outputs
    %   pred = guess of the position s occupies 
    %   in the L columns of V
    
   
    y = max(0,s'*W - Th); % compute reaction to s
    idx = find(y > 0); % find which neurons respond at all
    F = (V > Th)';
    numstim = sum(F,1);
    pred = 0;
    
    
    [~,l] = size(idx);
    for i=1:l
        if numstim(idx(i)) == 1 % find if one of the neurons is selective
            pred = find(F(:,idx(i)) == 1); % predict what the selective neuron says
            return;
        end
    end
    
    % if there are no neurons selective to s
    [~,best] = size(V);
    ind = 0;
    for i=1:l
        if numstim(idx(i)) < best % find the neuron that activates to less stimuli
            ind = i;
            best = numstim(idx(i));
        end
    end
    

    if ind ~= 0
        r = randi([1, best],1,1); % choose randomly from the stimuli it is active to
        fin = find(F(:,idx(i)) == 1,best);
        pred = fin(best);
        return;
    end
    
    % if there are no active neurons (i.e. s is a lost stimulus)
    % choose randomly from lost stimuli
    numneur = sum(F,2);
    idx = find(numneur == 0);
    [l,~] = size(idx);
    r = randi([1,l],1,1);
    pred = idx(r);
   
end 
    
    