%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sub Figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath functions\
%%

load TabSpace2.mat

figure
subplot(2,3,1)
imagesc(kkd1,kkb1,tabSkbn)
ylim([.001 1])
xlim([.5 10])
title("Parameter space")
xlabel('n')
ylabel('k_b')
set(gca,'ydir','normal','FontSize',12,'YScale','log','YMinorTick','on','XMinorTick','on')
pbaspect([1 1 1])
colormap(sky(3))
drawnow

subplot(2,3,2)
imagesc(kkd2,kkb2,tabSkbkd)
ylim([.001 1])
xlim([.01 1])
title("Parameter space")
xlabel('k_c')
ylabel('k_b')
set(gca,'ydir','normal','FontSize',12,'XScale','log','YScale','log')
colormap(sky(3))
pbaspect([1 1 1])
drawnow

subplot(2,3,3)
imagesc(kkd3,kkb3,tabSkcn)
ylim([.01 1])
xlim([.5 10])
title("Parameter space")
xlabel('n')
ylabel('k_c')
set(gca,'ydir','normal','FontSize',12,'YScale','log')
colormap(sky(3))
pbaspect([1 1 1])
drawnow

subplot(2,3,4)
imagesc(kkd4,kkb4,tabSKn)
ylim([.01 1])
xlim([.5 10])
title("Parameter space")
ylabel('K')
xlabel('n')
set(gca,'ydir','normal','FontSize',12,'YScale','log','XScale','linear')
colormap(sky(3))
pbaspect([1 1 1])
drawnow

subplot(2,3,5)
imagesc(kkd5,kkb5,tabSKkc)
ylim([.01 1])
xlim([.01 1])
title("Parameter space")
xlabel('K')
ylabel('k_c')
set(gca,'ydir','normal','FontSize',12,'XScale','log','YScale','log')
colormap(sky(3))
pbaspect([1 1 1])
drawnow

subplot(2,3,6)
imagesc(kkd6,kkb6,tabSKkb)
ylim([.001 1])
xlim([.01 1])
title("Parameter space")
xlabel('K')
ylabel('k_b')
set(gca,'ydir','normal','FontSize',12,'XScale','log','YScale','log')
colormap(sky(3))
pbaspect([1 1 1])
drawnow