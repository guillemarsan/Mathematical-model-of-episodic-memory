function pred = predictcon(W,U,s,Th,Thcn)
    % Function which, given a stimulus s, tries to predict
    % the class of it with the concept layer
    %
    % Inputs
    %   W = nxM weights matrix from stimulus to selective layer
    %   U = Mxa weights matrix from selective to concept layer
    %   s = nx1 stimulus shown
    %   Th = threshold selective layer
    %   Thcn = threshol concept layer
    % Outputs
    %   pred = guess of the class
    
    y = max(0,W'*s - Th);   % compute selective activation
    ycn = max(0, U'*y - Thcn);  % compute concept activation
    k = find(ycn);
    if isempty(k) || length(k) > 1
        pred = -1;   % return uncertainity
    else
        pred = k;   % return guess
    end
end


