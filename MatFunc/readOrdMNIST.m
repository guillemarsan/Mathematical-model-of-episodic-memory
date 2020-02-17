function [imgs, labls] = readOrdMNIST(num)
    % Function that outputs a 
    
    offset = 0;
    K = num/10;
    % read the database
    [d, l] = readMNIST("train-images.idx3-ubyte","train-labels.idx1-ubyte", num*5, offset);
    [~,idx] = sort(l); % sort just the first column
    d = d(:,:,idx);
      
    imgs = zeros(20,20,num);
    cont = 1;
    for i = 0:9
        imgs(:,:,i*K+1:(i+1)*K) = d(:,:,cont:cont+K-1);
        cont = cont + sum(l(:) == i);
        labls(i*K+1:(i+1)*K) = i;
    end   
     
end