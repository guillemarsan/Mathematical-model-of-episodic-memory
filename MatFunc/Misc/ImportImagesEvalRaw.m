function [imgs, class] = ImportImagesEvalRaw(FLDR, Figures, PlotFLG)

nIm = length(Figures);

imgs = [];
class = [];
if PlotFLG, cla; hold off; end

for k = 1:nIm
    file_pattern = fullfile(FLDR,[Figures{k}, '*.jpeg']);
    fls = dir(file_pattern);
    nFls = length(fls);
    for j = 1:nFls
        Img = imread(fullfile(FLDR, fls(j).name));                
        p = imbox(Img,0);
        
        imgs = [imgs, p(:)];

        if PlotFLG
            subplot(nIm,nFls, (k-1)*nFls + j)
            imagesc(Img); axis square; axis off
        end
    end
    class = [class, k*ones(1,length(fls))];
end

