function O = orderRasterPlot(F)
    s = sum(F,1);
    S = F(:,s == 1);
    A = F(:,s > 1);
    I = F(:,s == 0);

%     i = sum(S,2);
%     S = sortrows([i,S]);
%     S = S(:,2:end);

    i = sum(A,1);
    A = sortrows([i;A]')';
    A = A(2:end,:);

    [~,i] = max(S);
    S = sortrows([i;S]')';
    S = S(2:end,:);
    
    O = flip([S A I]);
   
end

