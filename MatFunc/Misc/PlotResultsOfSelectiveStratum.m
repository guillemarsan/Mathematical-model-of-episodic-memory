function PlotResultsOfSelectiveStratum(s, Resp0, Resp)

clf
subplot(3,1,1)
plot(s)
xlabel('dimension j','FontSize',18)
ylabel('s_j','FontSize',18)
xlim([0 size(s,1)])


subplot(3,2,3)
imagesc(Resp0)
colormap(flipud(colormap('gray')))
xlabel('neurons','FontSize',16)
ylabel('stimulus','FontSize',16)
title('neuronal response before learning','FontSize',18)

subplot(3,2,4)
imagesc(Resp)
colormap(flipud(colormap('gray')))
xlabel('neurons','FontSize',16)
ylabel('stimulus','FontSize',16)
title('neuronal response after learning','FontSize',18)

subplot(3,2,5)
R0 = sum(Resp0);
plot(R0)
xlim([0, size(Resp0,2)])
xlabel('neurons','FontSize',16)
ylabel('# responses','FontSize',16)
title(['Select. = ' num2str(sum(R0==1)) '; Inact. = ' num2str(sum(R0==0))...
    '; Rmean = ' num2str(mean(R0),2)]) 

subplot(3,2,6)
R = sum(Resp);
plot(R)
xlim([0, size(Resp,2)])
xlabel('neurons','FontSize',16)
ylabel('# responses','FontSize',16)
title(['Select. = ' num2str(sum(R==1)) '; Inact. = ' num2str(sum(R==0))...
     '; Rmean = ' num2str(mean(R),2)])