function ReImg = imbox2(Im, imageSize, plotop)
    % Function which finds the red element of the image and bounds it in a 
    % imageSize square
    %
    % Inputs
    %   Im = RGB matrix that represents the image
    %   imageSize = size of the output image
    %   plotop = option to plot the transformations
    % Outputs
    %   ReImg = RGB matrix that represents the cropped imageSize image
    
    % Color mask
    redness = double(Im(:,:,1)) - max(double(Im(:,:,2)),double(Im(:,:,3)));
    mask = redness > 50;
    BiImg = ones(size(Im,1),size(Im,2));
    BiImg = BiImg.*mask;
    FillBiImg = imfill(BiImg);

    % Bounding box
    sbox = regionprops(FillBiImg,'BoundingBox');
    box = floor(sbox.BoundingBox);  
    % Crop original pic
    BoxIm = Im(box(2)-30:box(2)+box(4)+30,box(1)-30:box(1)+box(3)+30,:);
    
    % Size required for the CNN
    ReImg = imresize(BoxIm,[imageSize(1) imageSize(2)]);
    if plotop
        figure();
        imshow(BiImg);
        figure();
        imshow(FillBiImg);
        figure();
        imshow(BoxIm);
        figure();
        imshow(ReImg);
    end    
end
