function [W,U] = SimulateNeurons4Layers(Tmax, h, W, U, x, f, K, a, b2, bcn2, Th, Thcn, d);

N = round(Tmax/h); % number of integration steps
NK = round(K/h); % number of integration steps that neurons integrate
Nsec = round(1/h); % number of integration steps of 1 sec

dS = zeros(size(W,2),1); % set new term to 0 at t = 0
dScn = zeros(size(U,2),1);
scon = 0;
scon2 = 0;
for k = 1:N-1
    % Predictor selective
    ii = f(k*h);      % find active stimulus at t = kh    
    [dw, y] = EvRHS(W, x(:,ii), Th, b2, dS);
    p = W + h*a*dw;
    
    % Predictor concept
    if mod(k,NK) == 1
        scon = y;
    elseif mod(k,Nsec) == 1
        scon = scon+y;
    end
    [du,ycn] = EvRHS(U, scon, Thcn, bcn2, dScn);
    pcn = U + h*a*du;
    
    % Corrector selective
    ii = f((k+1)*h);  % active stimulus at t = (k+1)h
    [dw2, y2] = EvRHS(p, x(:,ii), Th, b2, dS);
    W = W + 0.5*h*a*(dw + dw2); 
    
    % Corrector concept
    if k == 1
        scon2 = scon;
    elseif mod((k+1),NK) == 1
        scon2 = y2;
    elseif mod((k+1),Nsec) == 1
        scon2 = scon2+y2;
    end
    [du2,~] = EvRHS(pcn, scon2, Thcn, bcn2, dScn);
    U = U +0.5*h*a*(du + du2);
    
    % Update new term selective
    vin = d*(sum(y.^2) - sum(y)*y);
    dS = max(0, vin - Th);
    
    % Update new term concept
    vincn = d*(sum(ycn.^2) - sum(ycn)*ycn);
    dScn = max(0, vincn - Thcn);
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











