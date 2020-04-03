function showPattern(p)
    % Function to show patterns
    % Inputs
    %   p = kxk binary matrix
    
    [r, c] = size(p);                          % Get the matrix size
    imagesc((1:c)+0.5, (1:r)+0.5, p);
    colormap(flipud(gray));                              % Use a gray colormap
    axis equal
    set(gca, 'XTick', 1:(c+1), 'YTick', 1:(r+1), ...  % Change some axes properties
         'XLim', [1 c+1], 'YLim', [1 r+1], ...
         'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');
end