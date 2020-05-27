function pred = predictcon(W,U,s,Th,Thcn,dict)
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
    vot = dict(ycn > 0);
    con = unique(vot);
    if con == -1
        pred = -1;   % return uncertainity
    else
        pred = mode(vot(vot ~= -1));   % return guess
    end
end


