function pred = predictcon4(W,U,s,Th,Thcn,dict)
    % Function which, given a stimulus s, tries to predict
    % the class of it with the concept layer
    %
    % Inputs
    %   W = nxM weights matrix from stimulus to selective layer
    %   U = Mxa weights matrix from selective to concept layer
    %   s = nx1 stimulus shown
    %   Th = threshold selective layer
    %   Thcn = threshol concept layer
    %   dict = 1xa array that tags each conceptual neuron to the concept
    % Outputs
    %   pred = guess of the class
    
    y = max(0,W'*s - Th);   % compute selective activation
    ycn = max(0, U'*y - Thcn);  % compute concept activation
    concepts = length(unique(dict));
    if concepts == 1
        pred = -1;
    else
        pond = zeros(1,concepts-1);
        %pond(1) = sum(ycn(dict == -1));
        for i = 1:concepts-1
            pond(i) = sum(ycn(dict == i)); 
        end
        if sum(pond) ~= 0
            [~,indx] = max(pond);
        else
            indx = -1;
        end
        pred = indx;   
    end 
end