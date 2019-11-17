function ac = accuracy(W,s,Th)
    % Function that computes the error rate
    %
    % Inputs
    %   W = weight matrix nxM
    %   s = stimulus matrix nxL
    %   Th = parameter theta
    % Outputs
    %   ac = hits/L
    
    
    [~,L] = size(s);
    hits = 0;
    for i =1:L
        if(predict(W'*s,W,s(:,i),Th) == i)
            hits = hits + 1;
        end
    end
    ac = hits/L;
end