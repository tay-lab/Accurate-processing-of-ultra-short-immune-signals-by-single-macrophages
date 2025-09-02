%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath functions\
load AnalysisTotal.mat
%% Figure 1c
gname1=gname([3,5,10,15,20]);
gnameShort=["CpG","LPS","Pam","R848","TNF"];
Locb=ismember(data0.Category,gname1);
data=data0(Locb,:);

tspan=[-12 240];
clims=[0 6];            %%limits for NFkB amplitud
pp=[2 5];
plot_Fig1c(data,gname1,gnameShort,tspan,clims,pp,cc)