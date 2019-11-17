function W = SimulateNeurons2(Tmax, h, W, x, f, a, b2, Th, inh)

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
    dS = inh(y);
end

end

function [dw, y] = EvRHS(w, x, Th, b2, dS)
x = sum(x,2);   % sum of several simultaneous inputs (if exist)
v = w'*x;       % membrane potential
y = max(0, v - Th - dS); % neuronal response

dw = zeros(size(w)); % derivative of weights (without alpha)
for k = 1:size(w,2)
    dw(:,k) = y(k)*(b2*x - v(k)*w(:,k));
end

if(sum(sum(isnan(dw))) > 0) %check is there are NaN in the weights
    disp('Error');
    return
end

end











