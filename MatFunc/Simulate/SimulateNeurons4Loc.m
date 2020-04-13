function W = SimulateNeurons4Loc(Tmax, h, W, x, f, a, b2, Th, d, loc)

K = round(Tmax/h); % number of integration steps

dS = zeros(size(W,2),1); % set new term to 0 at t = 0
for k = 1:K-1
    % Predictor
    ii = f(k*h);      % find active stimulus at t = kh    
    [dw, y] = EvRHS(W, x(:,ii), Th, b2, dS);
    p = W + h*a*dw;
    % Corrector
    ii = f((k+1)*h);  % active stimulus at t = (k+1)h   
    W = W + 0.5*h*a*(dw + EvRHS(p, x(:,ii), Th, b2, dS)); 
    % Update new term
    vin = zeros(size(y,1),1);

    % With hard boundary
%     for l = 1:loc:size(y,1)-loc+1
%       vin(l:(l+loc-1)) = d*(sum(y(l:(l+loc-1)).^2) - sum(y(l:(l+loc-1)))*y(l:(l+loc-1)));
%     end

    % With moving window
    vin = d*(movsum(y.^2,loc) - movsum(y,loc).*y);
   
    dS = max(0, vin - Th);
end

end

function [dw, y] = EvRHS(w, x, Th, b2, dS)
x = sum(x,2);   % sum of several simultaneous inputs (if exist)
v = w'*x;       % membrane potential
y = max(0, v - Th - dS); % neuronal response

% derivative of weights (without alpha)
dw = y' .* (b2*x - v'.*w);

if(sum(sum(isnan(dw))) > 0) %check is there are NaN in the weights
    disp('Error');
    return
end

end











