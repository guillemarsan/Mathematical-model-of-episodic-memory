function [e,mom1,ph,distances] = moments(Im, signa,plotop)
    % Function which computes the invariant moments of the red element 
    % of an image
    %
    % Inputs
    %   Im = RGB matrix that represents the image
    %   signa = number of points of radius signature to return
    %   plotop = option to plot the transformations
    % Outputs
    %   e = eccentricity
    %   mom1 = perimeter^2/area
    %   ph = Hu invaraints
    %   distances = array of signa elements of equally distributed points
    %   of the radius signature (with adaptation far-close)
    
    % Color mask
    redness = double(Im(:,:,1)) - max(double(Im(:,:,2)),double(Im(:,:,3)));
    mask = redness > 40;
    BiImg = ones(size(Im,1),size(Im,2));
    BiImg = BiImg.*mask;
    FillBiImg = imfill(BiImg);
    if plotop
        figure();
        imshow(BiImg);
        figure();
        imshow(FillBiImg);
    end

    % Eccentricity
    se = regionprops(FillBiImg,'Eccentricity');
    e = se.Eccentricity;

    % P^2/A
    p = regionprops(FillBiImg,'Perimeter');
    a = regionprops(FillBiImg,'Area');
    mom1 = p.Perimeter^2/a.Area;

    % Hu invariants
    ph = hu_invariants(FillBiImg);

    % Radius signature
    ce = regionprops(FillBiImg,'Centroid');
    c = ce.Centroid;
    boundary = bwboundaries(FillBiImg,'noholes');
    x = boundary{1}(:,2);
    y = boundary{1}(:,1);
    if plotop
        figure();
        imshow(FillBiImg);
        hold on
        scatter(c(1),c(2),'x','r');
        scatter(boundary{1}(:,2),boundary{1}(:,1),'.','r','LineWidth',0.001);
    end
    distances = sqrt((x - c(1)).^2 + (y - c(2)).^2);
    distances = distances / max(distances);
    
    % Pick only signa points
    pts = length(distances);
    step = 1/pts;
    distances = interp1(0:step:1-step,distances,0:1/signa:1-step);
    
    % Far-close adaptation
    distances(distances > 0.7) = 1;
    distances(distances < 0.7) = 0.5;
    distances(distances < 0.3) = 0;
    
    % Scale up
    distances = distances * 100;
end
