%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath functions\
load AnalysisTotal.mat

%% Figure 1c
gname1=gname([3,5,10,15,20]);
gname2=gname([2,4,9,14,19]);
gnameShort=["CpG","Lipid A","Pam","R848","TNF"];
Locb=ismember(data0.Category,gname1);
data=data0(Locb,:);
Locb=ismember(data0.Category,gname2);
data2=data0(Locb,:);

tspan=[-12 240];
clims=[0 6];            %%limits for NFkB amplitud
pp=[2 5];
plot_Fig1c(data,data2,gname1,gname2,gnameShort,tspan,clims,pp,cc)

% countlabels(data0.Category)
