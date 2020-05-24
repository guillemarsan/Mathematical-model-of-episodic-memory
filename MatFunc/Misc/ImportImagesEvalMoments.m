function [mom, class] = ImportImagesEvalMoments(FLDR, Figures, signa, PlotFLG)

nIm = length(Figures);


mom = [];
class = [];
if PlotFLG, cla; hold off; end

for k = 1:nIm
    file_pattern = fullfile(FLDR,[Figures{k}, '*.jpeg']);
    fls = dir(file_pattern);
    nFls = length(fls);
    for j = 1:nFls
        Img = imread(fullfile(FLDR, fls(j).name));                
        momi = moments(Img, signa, 0);
        
        mom = [mom, momi];

        if PlotFLG
            subplot(nIm,nFls, (k-1)*nFls + j)
            imagesc(Img); axis square; axis off
        end
    end
    class = [class, k*ones(1,length(fls))];
end

