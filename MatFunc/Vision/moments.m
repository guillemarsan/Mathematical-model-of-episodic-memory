function [e,mom1,ph,distances] = moments(Im, plotop)
% Im = imread(Im);
BiImg = ones(size(Im,1),size(Im,2));

% LAB
%rgb2lab
% ImF = im2double(Im);
% x = 4*ImF(:,:,1) - 2*ImF(:,:,2) - 2*ImF(:,:,3);
% [~,j] = max(max(x));
% [~,i] = max(x(:,j));
% maskr = Im(:,:,1) > Im(i,j,1) - 10; %90
% maskg = Im(:,:,2) < Im(i,j,1) + 5; %60
% maskb = Im(:,:,3) < Im(i,j,1) + 5; %60
% box = regionprops(FillBiImg,'BoundingBox')
% imshow(FillBiImg(512:861,701:1053))
% Ir = imrsize(ib,[20,20])
% Square
% if type == 0
%     maskr = Im(:,:,1) > 120;
%     maskg = Im(:,:,2) < 50;
%     maskb = Im(:,:,3) < 75;
% elseif type == 1 % Triangle
%     maskr = Im(:,:,1) > 90;
%     maskg = Im(:,:,2) < 60;
%     maskb = Im(:,:,3) < 60;
% elseif type > 1 % Rest
%     maskr = Im(:,:,1) > 160;
%     maskg = Im(:,:,2) < 80;
%     maskb = Im(:,:,3) < 80;
% end

% maskr = Im(:,:,1) > uint8(graythresh(Im(:,:,1))*255);
% maskg = Im(:,:,2) < uint8(graythresh(Im(:,:,2))*255);
% maskb = Im(:,:,3) < uint8(graythresh(Im(:,:,3))*255);
% mask = maskr .* maskg .* maskb;

% Color mask
redness = double(Im(:,:,1)) - max(double(Im(:,:,2)),double(Im(:,:,3)));
mask = redness > 50;
BiImg = BiImg.*mask;
FillBiImg = imfill(BiImg);
if plotop
    figure();
    imshow(BiImg);
    figure();
    imshow(FillBiImg);
end

% Biggest blob
% [labeledImage, ~] = bwlabel(FillBiImg,4);
% blobMeasurements = regionprops(labeledImage, 'area');
% allAreas = [blobMeasurements.Area];
% [~, sortIndexes] = sort(allAreas, 'descend');
% biggestBlob = ismember(labeledImage, sortIndexes(1));
% FillBiImg = biggestBlob > 0;
% if plotop
%     figure();
%     imshow(FillBiImg);
% end


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
[~,argmax] = max(distances);
pts = length(distances);
distances = circshift(distances,pts-argmax+1);
distances = distances / max(distances);
distances(distances > 0.5) = 1;
distances(distances < 0.5) = 0;
% distances = 1./(1 + exp(-1000.*(distances-0.5)));
% figure();
% scatter(1:pts,distances','.');
step = fix(pts/400);
offset = rem(pts,400);
distances = distances(1:step:pts-offset);
end
