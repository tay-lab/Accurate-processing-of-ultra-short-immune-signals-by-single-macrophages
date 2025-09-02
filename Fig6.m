%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath functions\
load FitParams.mat
%% Figure B
load Phasekb_n1.mat
col=sky(3);
figure
nn1=1.5;
kb1=[.00005,.0005,.005];
for i=3:-1:1
[my,~,brst]=getScore2(kb1(i),kd,nn1,kl,ka,0,1);

plot(brst,my,"LineWidth",3,"Color",col(i,:))
hold on
end
ylabel("Max Amplitude")
xlabel("Pulse Duration")
xticks(10.^[0:3])
xscale("log")
pbaspect([1 1 1])
set(gca,"FontSize",14)
legend(["m_1","m_2","m_3"])

%% Figure C
figure
load Phasekb_kd.mat
col=sky(3);
imagesc(kb,kd,b(:,:,1)',[-1 1])
colormap(col)
colorbar('YTick',linspace(-1,1,3))
ylim([.05 2])
yticks([0:.5:2])
xticks([0:.001:.005])
ylabel("k_d")
xlabel("k_b")
set(gca, 'YDir', 'normal',"FontSize",15);
text(.0001,1.8,"m_3","FontSize",15)
text(2e-3,1,"m_2","FontSize",15)
text(4e-3,.15,"m_1","FontSize",15)
pbaspect([1 1 1])
%% Figure D
figure
load Phasekb_n1.mat
imagesc(kb,n1,b(:,:,1)',[-1 1])
colormap(col)
colorbar('YTick',linspace(-1,1,3))
yticks([1:1:5])
xticks([0:.001:.005])
ylabel("n")
xlabel("k_b")
set(gca, 'YDir', 'normal',"FontSize",15);
text(1e-4,4,"m_3","FontSize",15)
text(2e-3,3,"m_2","FontSize",15)
text(4e-3,1,"m_1","FontSize",15)
pbaspect([1 1 1])
%% Figure E
figure
load Phasekb_K.mat
col=sky(3);
imagesc(kb,kl,b(:,:,1)',[-1 1])
colormap(col)
colorbar('YTick',linspace(-1,1,3))
yticks([0.00001:.00003:.0001])
xticks([0:.001:.005])
ylabel("K")
xlabel("k_b")
set(gca, 'YDir', 'normal',"FontSize",15);
text(1e-4,8e-5,"m_3","FontSize",15)
text(2e-3,5e-5,"m_2","FontSize",15)
text(4e-3,1e-5,"m_1","FontSize",15)
pbaspect([1 1 1])

%% Figure F
figure
load Phasekb_ka.mat
col=sky(3);
imagesc(kb,ka,b(:,:,1)',[-1 1])
colormap(col)
colorbar('YTick',linspace(-1,1,3))
yticks([0:.1:.5])
xticks([0:.001:.005])
ylabel("k_a")
xlabel("k_b")
set(gca, 'YDir', 'normal',"FontSize",15);
text(1e-4,1e-1,"m_3","FontSize",15)
text(2e-3,2e-1,"m_2","FontSize",15)
text(4e-3,4e-1,"m_1","FontSize",15)
pbaspect([1 1 1])
%% Figure G
load Phasekb_n1.mat
col=sky(3);
nn1=1.5;
kb1=[.00005,.0005,.005];
figure
for i=3:-1:1
    subplot(1,3,i)
[my,Mod,brst]=getScore2(kb1(i),kd,nn1,kl,ka,0,0);
for j=1:4
T=timetable(minutes(Mod{j}.Time),Mod{j}.YY(:,4));
T=synchronize(T,minutes([-10:200]),"linear");

plot(minutes(T.Time),T.Var1,"LineWidth",3)
hold on
end
ylabel("Nuclear NFkB")
xlabel("Pulse Duration")
set(gca,"FontSize",8)
ylim([0 1])
xlim([-10 180])
xticks([0:60:180])
pbaspect([1 1 1])
end
legend(["1","10","100","1000"])
