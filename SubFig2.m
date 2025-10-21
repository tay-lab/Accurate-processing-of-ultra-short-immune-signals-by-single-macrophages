%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SubFigure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
addpath functions\
load AnalysisTotal.mat
gname1=gname([3,5,10,15,20]);
gnameShort=["CpG","Lipid A","Pam","R848","TNF"];
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
%% Fig 2
%%Spatio-temporal Analysis
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
data.clsT=clsT;

%%
tt=data.Time(1,:);
kmax=3;
clim=[0 6];
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
end
end
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
if sp==5
colormap sky
   figure
   sp=1;
end
subplot(2,2,sp)
[rrM,pval]=corr(table2array(ERtable),table2array(LRtable));
rrM=rrM.^2;
rrM=round(rrM,2);
h=heatmap(rrM,"ColorLimits",[0 1]);
colormap turbo
h.XDisplayLabels=featLabel;
h.YDisplayLabels=featLabel;
xlabel("Early Responders")
ylabel("Late Responders")
title(gnameShort{j})
end


colormap sky