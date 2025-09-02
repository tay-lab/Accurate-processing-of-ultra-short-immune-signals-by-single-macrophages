%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
addpath functions\
load AnalysisTotal.mat
gname1=gname([3,5,10,15,20]);
gnameShort=["CpG","LPS","Pam","R848","TNF"];
Locb=ismember(data0.Category,gname1);
data=data0(Locb,:);

features=["EAUC","LAUC",...
        "Peak",...
        "Speed",'Fourier',...
        'Duration',"Activation","Time to Peak"];

featLabel=[" Early AUC"," Late AUC",...
        " Amplitude",...
        " Speed",' Oscillations',...
        ' Duration'," Activation Time", " Time to Peak"];

zz=data.Category==gname1{1};
data2=data(~zz,:);
datvio={};
datmn=[];
stdmn=[];
maxvio=[];
for ss=1:length(gname1)
    zz=data.Category==gname1(ss);
    for i=1:length(features)
        dd=data.(features(i))(zz);
        dd(dd<0)=0;
        datmn(ss,i)=median(dd,'omitnan');
        datvio{ss,i}=dd;
        mxvio(ss,i)=prctile(dd,99);
        stdmn(ss,i)=std(dd,'omitnan');
    end
end
col=hsv(6);
figure("Name","Violin Plots of Features Distributions")
kk=0;
for i=1:length(features)
    kk=kk+1;
    if i==10
        figure
        kk=1;
    end
    subplot(3,3,kk)
    dd=data2.(features(i));
    ss=data2.Category;
    zz2=isnan(dd);
    dd(zz2)=[];
    ss(zz2)=[];
    dd(dd<0)=0;
    violin({datvio{1,i},datvio{2,i},datvio{3,i},datvio{4,i},datvio{5,i}},...
        'facecolor',col,'mc',[],'medc','black','bw',mean(stdmn(:,i))/3)
    xticks(1:6)
    ylim([0 prctile(dd,99)])
    xticklabels(gnameShort(1:end))
    ylabel(featLabel(i))
    set(gca,'FontSize',10,'FontName','Times New Roman')
end

%% Fig 2
%%Spatio-temporal Analysis

figure
col=hsv(3);
kmax=3;
rng default

data2=data;

Nfkb=data2.NfkB;

ratings=[data2.Activation data2.("Peak")];
ratings2=ratings;

nzn=isnan(ratings(:,1));
nzn2=isnan(ratings(:,2));
ratings(nzn | nzn2,1)=1000;
ratings(nzn | nzn2,2)=1000;
ratings2(nzn | nzn2,:)=[];

[clsT,cnt] = kmeans(ratings,kmax,'distance','cityblock','MaxIter',100,'Replicates',10);

auc=[];
for i=1:kmax
    auc=[auc; mean(ratings(clsT==i,1))];
end
[~,idx]=sort(auc,"ascend");

clsTp=clsT;
for i=1:kmax
    clsT(clsTp==idx(i))=i;
end
gscatter(ratings(:,1),ratings(:,2),clsT,col,'...',[16 16 16])
data.clsT=clsT;

xlim([0 250])
ylim([0 20])

xlabel("Activation time")
ylabel("Peak")
set(gca,'FontSize',14)
legend({"Early Responders","Late Responders"})

%%
tt=data.Time(1,:);
kmax=3;
clim=[0 6];
figure
sps=0;
pCl=zeros(5,kmax);
rgname=[" First Responders"," Second Responders"," Non-Responders"];
for j=1:5
    sps=sps+1;
zz=data.Category==gname1{j};
data2=data(zz,:);
NfkB=data2.NfkB;

clsT=data2.clsT;
    nm=[];
for i=1:kmax
    pl=nnz(clsT==i)/nnz(clsT);
    pl=round(100*pl,0);
    pCl(j,i)=pl;
    nnf=NfkB(clsT==i,:);
    auc=data2.Activation(clsT==i,:);
    [auc,la]=sort(auc,"ascend");

    if j==4 & i==1
        sps=1;
        figure
    end
        subplot(3,kmax,i+(sps-1)*kmax)

    imagesc(tt,1:size(nnf,1),nnf(la,:),clim)
    colormap parula


nm=["f="+num2str(pl)+"%"];
xlim([-10 180])
xticks([0:60:240])
xlabel("Time (min)")
ylabel([])
yticks([])
title(gnameShort{j}+rgname(i)+" "+nm,"FontSize",15)
set(gca,'FontSize',10,'FontName','Times New Roman')
end
end

%%
figure
bar(gnameShort,pCl)
legend({"ER","LR","NR"})
set(gca,'FontSize',14)
%%

figure
rr=[];
sp=0;
for j=1:5
    sp=sp+1;
zz=data.Category==gname1{j};
data2=data(zz,:);
gr=unique(data2.cellNum(:,1));
mmd1=[];
mmd2=[];
mmd3=[];
mmd4=[];
acTime=[];
acLAUC=[];
acPeak=[];
ERtable=[];
LRtable=[];
NFtable=[];

for i=1:length(gr)
    zz=string(data2.cellNum(:,1))==gr{i};
    data3=data2(zz,features);
    tempX=data2.X(zz,:);
    tempY=data2.Y(zz,:);
    Nf=data2.NfkB(zz,:);

    aTime=data2.Activation(zz);
    aEAUC=data2.EAUC(zz);
    aPeak=data2.Peak(zz);
    aLAUC=data2.LAUC(zz);

    cT=data2.clsT(zz);
    fr=find(cT==1);
    sr=find(cT==2);

    if ~isempty(fr) && ~isempty(sr)
    TotTime=find(tt>=0 & tt<=45);
    srC=zeros(numel(sr),numel(TotTime),4);
    frC=zeros(numel(sr),numel(TotTime));

    for tim=TotTime
        XY=[tempX(:,tim) tempY(:,tim)];
        [Idx,Ds] = knnsearch(XY(fr,:),XY(sr,:));
        srC(:,tim-TotTime(1)+1,1)=Ds;
        frC(:,tim-TotTime(1)+1)=fr(Idx);

        srC(:,tim-TotTime(1)+1,2)=aEAUC(fr(Idx));
        srC(:,tim-TotTime(1)+1,3)=aPeak(fr(Idx));
        
        [Idx,Ds] = rangesearch(XY(fr,:),XY(sr,:),300);
        ie=cellfun(@length,Ds);
        srC(:,tim-TotTime(1)+1,4)=ie;
    end
    [~,idx]=min(srC(:,:,1),[],2,"linear");

    ERtable=[ERtable; data3(frC(idx),:)];
    LRtable=[LRtable; data3(sr,:)];
    NFtable=[NFtable; Nf(frC(idx),:)];

    mmd1=[mmd1; mean(srC(:,:,1),2)];
    mmd2=[mmd2; mean(srC(:,:,2),2)];
    mmd3=[mmd3; mean(srC(:,:,3),2)];
    mmd4=[mmd4; median(srC(:,:,4),2)];

    acTime=[acTime; aTime(sr,:)];
    acLAUC=[acLAUC; aLAUC(sr,:)];

    end
    
end
% 
if sp==6
   figure
   sp=1;
end


subplot(4,5,sp)
scatter(mmd1,acTime)
rr=corr(mmd1,acTime);
rr=rr.^2;

text(prctile(mmd1,80),prctile(acTime,95),"r^2="+num2str(round(rr,2)),'FontSize',14)
l=lsline;
set(l,"LineWidth", 2,"Color","r")

title(gnameShort{j})
% xlim([min(mmd1) max(mmd1)+5])
% ylim([min(acTime)-5 max(acTime)+5])
xlabel("Distance to ER")
ylabel("Activation Time")
set(gca,'FontSize',8)

subplot(4,5,sp+5)
scatter(mmd2,acLAUC)
rr=corr(mmd2,acLAUC);
rr=rr.^2;

l=lsline;
set(l,"LineWidth", 2,"Color","r")
text(prctile(mmd2,80),prctile(acLAUC,95),"r^2="+num2str(round(rr,2)),'FontSize',14)
% xlim([min(mmd2) max(mmd2)+5])
% ylim([min(acTime)-5 max(acTime)+5])
xlabel("EAUC")
ylabel("LAUC")
set(gca,'FontSize',8)

subplot(4,5,sp+10)
scatter(mmd3,acLAUC)
rr=corr(mmd3,acLAUC);
rr=rr.^2;

l=lsline;
set(l,"LineWidth", 2,"Color","r")
text(prctile(mmd3,80),prctile(acLAUC,95),"r^2="+num2str(round(rr,2)),'FontSize',14)
% xlim([min(mmd2) max(mmd2)+5])
% ylim([min(acTime)-5 max(acTime)+5])
xlabel("Peak")
ylabel("LAUC")
set(gca,'FontSize',8)

subplot(4,5,sp+15)
scatter(mmd4,acTime)
rr=corr(mmd4,acTime);
rr=rr.^2;

l=lsline;
set(l,"LineWidth", 2,"Color","r")
text(prctile(mmd4,80),prctile(acTime,95),"r^2="+num2str(round(rr,2)),'FontSize',14)
% xlim([min(mmd2) max(mmd2)+5])
% ylim([min(acTime)-5 max(acTime)+5])
ylabel("Activation Time")
xlabel("Num Neighbors")
set(gca,'FontSize',8)

end