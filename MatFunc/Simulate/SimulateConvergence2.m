function y2 = SimulateConvergence2(W,s,Th,d,eps,maxiter)
    

    V = W'*s;
    y = ones(size(W,2),size(s,2));
    y2 = max(0, V - Th);
    j = 0;
    while (norm(y - y2) > eps && j < maxiter)
        y = y2;
        vin = d*(sum(y.^2) - sum(y).*y);
        dS = max(0, vin - Th);
        y2 = max(0, V - Th - dS);
        j = j+1;
    end
    
    if j == maxiter
        print('No convergence');
    end
end

