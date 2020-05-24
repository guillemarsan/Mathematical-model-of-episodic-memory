function dict = conceptmap(F,K)
    % Function which, given an firing raster, and an integration step K
    % computes an array tagging conceptual neurons
    %
    % Inputs
    %   F = axL firing raster last layer neurons - stimuli
    %   K = integration step
    % Outputs
    %   dict = 1xa array where each value is the concept each conceptual
    %   neuron associates. -1 if it is not conceptul
    
    [A,L] = size(F);
    dict = ones(1,A)*(-1);
    % For each concept
    for i = 0:(L/K)-1
        onesidx = i*K+1:(i+1)*K;
        zerosidx = setdiff(1:L,onesidx);
        zerarr = sum(F(:,zerosidx),2);
        onearr = sum(F(:,onesidx),2);
        % If it activates only for all stimuli of the concept
        dict((zerarr == 0) & (onearr == K)) = i;
    end

end

