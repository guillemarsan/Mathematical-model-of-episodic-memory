function W = SimulateNeurons3(Tmax, h, W, x, f, a, b2, Th, d)

K = round(Tmax/h); % number of integration steps
for k = 1:K-1
    % Predictor
    ii = f(k*h);      % find active stimulus at t = kh
    dw = EvRHS(W, x(:,ii), Th, b2, d);
    p = W + h*a*dw;
    % Corrector
    ii = f((k+1)*h);  % active stimulus at t = (k+1)h
    W = W + 0.5*h*a*(dw + EvRHS(p, x(:,ii), Th, b2, d)); 
end

end

function dw = EvRHS(w, x, Th, b2, d)
x = sum(x,2);   % sum of several simultaneous inputs (if exist)
v = w'*x;       % membrane potential
y = max(0, v - Th); % neuronal response
dw = zeros(size(w)); % derivative of weights (without alpha)

dum = b2 - d*(sum(y) - y); % new term
for k = 1:size(w,2)
    dw(:,k) = y(k)*(dum(k)*x - v(k)*w(:,k));
end
end











