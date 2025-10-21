%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath functions\
load AnalysisTotal.mat
%%
gname1=gname([1,3,5,10,15,20]);
gnameShort=["Control","CpG","Lipid A","Pam","R848","TNF"];
Locb=ismember(data0.Category,gname1);
data=data0(Locb,:);
%%
MIm=zeros(6);
TTm=zeros(6);
tTm=zeros(6);
    tt=[data.("Time")(1,:) ];
for kk=1:6
for j=1:6
    zz1=data.Category==gname1(kk);
    zz2=data.Category==gname1(j);

    S=data.Category;
    Nfkb=[data.("NfkB")];
    k=10;
    MIt=zeros(size(Nfkb,2),1);
    parfor i=1:size(Nfkb,2)
        Rsp=[Nfkb(zz1,i); Nfkb(zz2,i)];
        Sp=[S(zz1); S(zz2)];
        MIt(i)=getCC(Rsp,Sp,k);
    end
    [mm,tm]=max(MIt);
    MIm(kk,j)=mm;
    TTm(kk,j)=tt(tm);
    tTm(kk,j)=tm;
end
end

%%
figure
colormap sky
MIm=round(triu(MIm,1),2);
MIm(MIm<=0)=NaN;
h=heatmap(MIm, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
h.XDisplayLabels=string(gnameShort);
h.YDisplayLabels=string(gnameShort);
h.ColorLimits = [0 1];
set(gca,'XLim',[2,6],'yLim',[1,5])
set(gca,'FontSize',18,'FontName','Times New Roman')
title("Ligand specificity")

figure
TTm=round(triu(TTm,1),2);
TTm(TTm<=0)=NaN;
h=heatmap(TTm, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
h.XDisplayLabels=string(gnameShort);
h.YDisplayLabels=string(gnameShort);
h.ColorLimits = [0 70];
set(gca,'XLim',[2,6],'yLim',[1,5])
set(gca,'FontSize',18,'FontName','Times New Roman')
title("Time to Max specificity")
colormap summer

%%
S=data.Category;
Nfkb=[data.("NfkB")];
tt=[data.("Time")(1,:)];
k=10;
MI=zeros(size(Nfkb,2),1);

parfor i=1:size(Nfkb,2)
    Rs=Nfkb(:,i);
    MI(i)=getCC(Rs,S,k);
end

figure("Name","Channel capacity")
[pk,lpk]=findpeaks(MI,"SortStr","descend","NPeaks",1);

plot(tt,MI,'LineWidth',2)
ylabel("Channel Capacity (bits)")
xlabel("Time (min)")
xlim([-10 240])
ylim([0 1])
hold on
plot([tt(lpk),tt(lpk)], [0,pk],'LineStyle','--','LineWidth',2)
xticks([0,tt(lpk),60:60:240])
xticklabels([0,"t_{max}",60:60:180])
text(tt(lpk),pk,"NFkB CC_{max}",'HorizontalAlignment','left',...
    'VerticalAlignment','bottom','FontSize',14,'FontName','Times New Roman')
set(gca,'FontSize',14,'FontName','Times New Roman')
hold off
title("Channel Capacity")

%%
data.("CC NFkB")=Nfkb(:,lpk);
features=["CC NFkB","EAUC","LAUC",...
        "Peak",...
        "Speed",'Fourier',...
        'Duration',"Activation","Time to Peak"];

MIfeat=zeros(length(features),1);
parfor i=1:length(features)
    Rs=data.(features(i));
    MIfeat(i)=getCC(Rs,S,k);
end

figure("Name","Feature score")
featLabel=[" NFkB CC_{max}"," Early AUC"," Late AUC",...
        " Amplitude",...
        " Speed",' Oscillations',...
        ' Duration'," Activation", " Time to Peak"];


[~,idx1]=sort(MIfeat,"descend");

bh=barh(MIfeat(idx1(1:end)),'r');
set(gca,'YDir','reverse','Xlim',[0 1.4])
set(gca,'FontSize',14,'FontName','Times New Roman')
yticks([])
ylabel("Features")
xlabel("Mutual Information (bits)")
title("Information storaged in NFkB features")

ytips1 = bh(1).XEndPoints;
xtips1 = bh(1).YEndPoints;
text(xtips1,ytips1,featLabel(idx1(1:end)),'HorizontalAlignment','left',...
    'VerticalAlignment','middle','FontSize',14,'FontName','Times New Roman')

%%
figure("Name","Ligand specificity")
nlas=[1,2,4];
MIcomb=[];
k=10;
Nfkb=[data.("NfkB")];
for i=1:length(unique(S))
    zz1=data.Category==gname1(i);
    for j=1:length(unique(S))
        zz2=data.Category==gname1(j);
        Rs=[Nfkb(zz1,tTm(i,j)) table2array(data(zz1,features(idx1(nlas)))); ...
            Nfkb(zz2,tTm(i,j)) table2array(data(zz2,features(idx1(nlas))))];
        nc1=sum(zz1);
        nc2=sum(zz2);
        Ss=categorical([ones(nc1,1); 2*ones(nc2,1)]);
        MIcomb(i,j)=getCC(Rs,Ss,k);
    end
end

MIcomb=round(triu(MIcomb,1),2);
MIcomb(MIcomb<=0)=NaN;
h=heatmap(MIcomb, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
h.XDisplayLabels=string(gnameShort);
h.YDisplayLabels=string(gnameShort);
h.ColorLimits = [0 1];
set(gca,'XLim',[2,6],'yLim',[1,5])
set(gca,'FontSize',18)
title("Ligand specificity")
colormap sky
%%
% save functions/AnalysisTotal data0 data groups gname cc tspan sm ActiveTresh cnt MIm TTm tt