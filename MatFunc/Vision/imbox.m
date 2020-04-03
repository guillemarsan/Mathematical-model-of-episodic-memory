function ReBFBImg = imbox(Im, plotop)
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
% if type == 0
%     maskr = Im(:,:,1) > 120;
%     maskg = Im(:,:,2) < 50;
%     maskb = Im(:,:,3) < 75;
% elseif type == 1
%     maskr = Im(:,:,1) > 90;
%     maskg = Im(:,:,2) < 60;
%     maskb = Im(:,:,3) < 60;
% elseif type > 1
%     maskr = Im(:,:,1) > 140;
%     maskg = Im(:,:,2) < 80;
%     maskb = Im(:,:,3) < 80;
% end
% mask = maskr .* maskg .* maskb;

% Color mask
redness = double(Im(:,:,1)) - max(double(Im(:,:,2)),double(Im(:,:,3)));
mask = redness > 50;
BiImg = BiImg.*mask;
FillBiImg = imfill(BiImg);

% Bounding box
sbox = regionprops(FillBiImg,'BoundingBox');
box = floor(sbox.BoundingBox);
BoxFBImg = FillBiImg(box(2)-20:box(2)+box(4)+20,box(1)-20:box(1)+box(3)+20);
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
