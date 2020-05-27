function pred = predictcon2(W,U,s,Th,Thcn,d,dict)
    % Function which, given a stimulus s, tries to predict
    % the class of it with the concept layer. Uses inhibition convergence 
    % to compute activation
    %
    % Inputs
    %   W = nxM weights matrix from stimulus to selective layer
    %   U = Mxa weights matrix from selective to concept layer
    %   s = nx1 stimulus shown
    %   Th = threshold selective layer
    %   Thcn = threshold concept layer
    %   d = inhbition parameter
    %   dict = 1xa array that tags each conceptual neuron to the concept
    % Outputs
    %   pred = guess of the class
    
    y = SimulateConvergence2(W,s,Th,d,0.01,1000);
    ycn = SimulateConvergence2(U,y,Thcn,d,0.01,1000);
    vot = dict(ycn > 0);
    con = unique(vot);
    if con == -1
        pred = -1;   % return uncertainity
    else
        pred = mode(vot(vot ~= -1));   % return guess
    end
end


