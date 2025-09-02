%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SubFigure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load AnalysisTotal.mat
addpath functions\
%%
gname1=gname([1,2,4,9,14,19]);
gnameShort=["Control","CpG","LPS","Pam","R848","TNF"];
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
ylim([0 1.2])
hold on
plot([tt(lpk),tt(lpk)], [0,pk],'LineStyle','--','LineWidth',2)
xticks([0,60,tt(lpk),120:60:240])
xticklabels([0,"t_{max}",60:60:180])
text(tt(lpk),pk,"NFkB CC_{max}",'HorizontalAlignment','left',...
    'VerticalAlignment','bottom','FontSize',14,'FontName','Times New Roman')
set(gca,'FontSize',14,'FontName','Times New Roman')
hold off
title("Channel Capacity")