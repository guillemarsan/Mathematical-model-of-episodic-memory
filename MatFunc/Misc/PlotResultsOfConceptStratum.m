function PlotResultsOfConceptStratum(y, Resp0, Resp)

clf
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
subplot(3,1,1)
plot(y)
xlabel('dimension j','FontSize',18)
ylabel('y_j','FontSize',18)
xlim([0 size(y,1)])


subplot(3,2,3)
R0 = orderRasterPlot(Resp0);
imagesc(R0)
colormap(flipud(colormap('gray')))
xlabel('neurons','FontSize',16)
ylabel('stimulus','FontSize',16)
title('neuronal response before learning','FontSize',18)

subplot(3,2,4)
R = orderRasterPlot(Resp);
imagesc(R)
colormap(flipud(colormap('gray')))
xlabel('neurons','FontSize',16)
ylabel('stimulus','FontSize',16)
title('neuronal response after learning','FontSize',18)

dict0 = conceptmap(Resp0',3);
subplot(3,2,5)
histogram(dict0,-1.25:0.5:3.25)
xlabel('types of neurons','FontSize',16)
ylabel('# neurons','FontSize',16)
title(['ConScore3 = ' num2str(conceptscore(dict0,3))...
    ';Concept 1 = ' num2str(sum(dict0==1)) '; Concept 2 = ' num2str(sum(dict0==2))...
    '; Concept 3 = ' num2str(sum(dict0==3)) '; Nonconceptual(-1) = ' num2str(sum(dict0==-1))]) 

dict = conceptmap(Resp',3);
subplot(3,2,6)
histogram(dict,-1.25:0.5:3.25)
xlabel('types of neurons','FontSize',16)
ylabel('# neurons','FontSize',16)
title(['ConScore3 = ' num2str(conceptscore(dict,3))...
    '; Concept 1 = ' num2str(sum(dict==1)) '; Concept 2 = ' num2str(sum(dict==2))...
    '; Concept 3 = ' num2str(sum(dict==3)) '; Nonconceptual(-1) = ' num2str(sum(dict==-1))]) 