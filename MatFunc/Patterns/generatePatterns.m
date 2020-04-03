function [p,s] = generatePatterns(k)
    % Function that generates simple patterns 
    % (lines) on a grid kxk
    %
    % Inputs
    %   k = dimension of the square grid
    % Outputs
    %   p = kxkx6k-2 matrix with the 6k-2 patterns
    %   s = k^2x6k-2 matrix. Is p with the grid 
    %   represented as a vector
    
p = zeros(k,k,6*k-2);
s = zeros(k*k,6*k-2);

% Horizontal and vertical lines
for i=1:k
    p(i,:,i) = ones(1,k);
    p(:,:,i+k) = p(:,:,i)';
end 

% Main diagonals
p(:,:,2*k+1) = diag(ones(1,k),0);
p(:,:,2*k+2) = flip(p(:,:,2*k+1));

% Smaller diagonals
for i=1:k-1
    p(:,:,i+2*k+2) = diag(ones(1,k-i),i);
    p(:,:,i+3*k+1) = p(:,:,i+2*k+2)'; % 2*k+2 + (k-1) = 3*k+1
    p(:,:,i+4*k) = flip(p(:,:,i+2*k+2));
    p(:,:,i+5*k-1) = flip(p(:,:,i+2*k+2),2);
end

% Vectorize the grids and normalize p
for i=1:6*k-2
    p(:,:,i) = p(:,:,i)/norm(p(:,:,i));
    aux = p(:,:,i)';
    s(:,i) = aux(:);
end
end


    