%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SubFigure 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath functions\
load AnalysisTotal.mat
%% Figure A
T=readtable('Plot Values.csv');
ti=T.Slice*.2-2;
tx=-2:.05:4;
signal=rescale(T.Mean);
gaussModel = 'a*exp(-((x-b)/c)^2/2)+d';
startPoints = [1.2 1 1 1];
[fittedModel, gof] = fit(ti, signal, gaussModel,'Start', startPoints);

cn = coeffnames(fittedModel);
cv = coeffvalues(fittedModel);
mn=cv(2);
sig=cv(3);

xconf=mn-sig:.01:mn+sig;
tsig=xconf(end)-xconf(1);

ti=ti-mn;

area(xconf-mn,fittedModel(xconf),"FaceColor",[1 .9 .9])
hold on
plot(tx-mn,fittedModel(tx),"black","LineWidth",2)
plot(ti,signal,"*r","LineWidth",2)
xlabel("Time (s)")
ylabel("Norm. Intensity")
xlim([-2 2])
xticks([-2:.5:2])
xticklabels([-1:.5:4])
legend(["65%","Data","Fit"])
%% Figure B
gname1=gname([1,2,4,9,14,19]);
gnameShort=["Control","CpG","LPS","Pam","R848","TNF"];
Locb=ismember(data0.Category,gname1);
data=data0(Locb,:);

tspan=[-12 240];
clims=[0 12];            %%limits for NFkB amplitud
pp=[2 6];
plot_Fig1c(data,gname1,gnameShort,tspan,clims,pp,cc)