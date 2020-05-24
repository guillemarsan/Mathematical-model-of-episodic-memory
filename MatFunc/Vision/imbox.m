function ReBFBImg = imbox(Im, plotop)
    % Function which finds the red element of the image and bounds it in a 
    % 20x20 square
    %
    % Inputs
    %   Im = RGB matrix that represents the image
    %   plotop = option to plot the transformations
    % Outputs
    %   ReBFBImg = 20x20 binary image that represents the red element of
    %   the image
    
    % Color mask
    redness = double(Im(:,:,1)) - max(double(Im(:,:,2)),double(Im(:,:,3)));
    mask = redness > 50;
    BiImg = ones(size(Im,1),size(Im,2));
    BiImg = BiImg.*mask;
    FillBiImg = imfill(BiImg);

    % Bounding box
    sbox = regionprops(FillBiImg,'BoundingBox');
    box = floor(sbox.BoundingBox);
    BoxFBImg = FillBiImg(box(2)-30:box(2)+box(4)+30,box(1)-30:box(1)+box(3)+30);
    ReBFBImg = imresize(BoxFBImg,[20,20]);
    if plotop
        figure();
        imshow(BiImg);
        figure();
        imshow(FillBiImg);
        figure();
        imshow(BoxFBImg);
        figure();
        showPattern(ReBFBImg);
    end
end
