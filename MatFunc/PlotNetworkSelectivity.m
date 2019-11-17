function PlotNetworkSelectivity(W0, W, s, Th, d)

R0   = s'*W0 > Th+1e-7; % response matrix at t = 0
Rend = s'*W > Th+1e-7;  % response matrix at t = Tmax
Rslctv = 100*mean(sum(Rend,1) == 1);      % ratio of selective neurons
Rlost = 100*mean(sum(Rend,2) == 0);       % ratio of lost stimuli
nLost =  sum(sum(Rend,2) == 0); % number of lost stimuli
nSel =  sum(sum(Rend,1) == 1);  % number of selective neurons
nSil =  sum(sum(Rend,1) == 0);  % number of silent neurons

disp('*************')
disp(['rSelective = ' num2str(Rslctv,3) '%; rLost = ' num2str(Rlost,3) '%'])


subplot(2,1,1)
nDet = [sum(R0); sum(Rend)]';  
nDet(nDet == 0) = 0.1;
h1 = semilogy(nDet, 'LineWidth', 2);
ylim([0.1 max(nDet(:))])
h1(1).Color = 0.7*[1 1 1];
ytk = get(gca,'YTick');
ytl = get(gca,'YTickLabel');
ytl(ytk == 0.1) = {'0'};
ytl(ytk == 1) = {'1'};
set(gca,'YTickLabel',ytl)
% grid on
xlabel('neuron','FontSize',14); 
ylabel('# detected stimuli','FontSize',14)
legend({'t = 0','t = T_{max}'},'FontSize',12,'location','best')
title(['Selective = ' num2str(nSel) '; Silent = ' num2str(nSil)], 'FontSize',12)

subplot(2,1,2)
nExc = [sum(R0,2) sum(Rend,2)];
h2 = plot(nExc, 'LineWidth', 2);
h2(1).Color = 0.7*[1 1 1];
xlabel('stimulus','FontSize',14); 
ylabel('# excited neurons','FontSize',14)
ylim([0 35])
title(['d = ' num2str(d) '; Lost = ' num2str(nLost)], 'FontSize',12)

