function y2 = SimulateConvergence1(W,s,Th,d,eps,maxiter, loc, hard)
    
    
    V = W'*s;
    y = ones(size(W,2),size(s,2));
    y2 = max(0, V - Th);
    j = 0;
    vin = zeros(size(V));
    while (norm(y - y2) > eps && j < maxiter)
        y = y2;
        if hard
            for l = 1:loc:size(y,1)-loc+1
                vin(l:(l+loc-1),:) = d*(sum(y(l:(l+loc-1),:).^2) - ...
                    sum(y(l:(l+loc-1),:)).*y(l:(l+loc-1),:));
            end
        else
            vin = d*(movsum(y.^2,loc) - movsum(y,loc).*y);
        end
        dS = max(0, vin - Th);
        y2 = max(0, V - Th - dS);
        j = j+1;
    end
    
    if j == maxiter
        disp('No convergence');
    end
end

