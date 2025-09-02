function plot_Fig1c(data,gname,gnameShort,tl,clims,p,cc)
rng(14)

figure("Name","Fig1")
count=0;


for chnum=1:numel(gname)
    if mod(count+1,17)==0
        h = axes('visible','off'); 
        c = colorbar(h,'Position',[0.93 0.12 0.022 0.7]);  % attach colorbar to h
        caxis(h,clims); 
        figure
        count=0;
    end
count=count+1;

dd=[data.("NfkB")];
ftime=[data.("Time")]; 
zz=data.Category==gname{chnum};
% zz=ismember(data.cellNum(:,1),gname{chnum});
ftime=ftime(zz,:);
dd=dd(zz,:);
% dd=smoothdata(dd',"gaussian",[3 3]);
ncells=size(dd,1);


% auc=trapz(dd,2);
auc=trapz(dd(:,1:30),2);
[~,Ind]=sort(auc,'descend');

    subplot(p(1),p(2),count)
    imagesc(ftime(1,:),1:ncells,dd(Ind,:),clims)
    xticks(0:120:tl(2))
    xlabel('Time (min)')
    yticks([])
    yticklabels([])
    ylim([0 ncells]);
    xlim([tl(1)  tl(2)]);
    title(gnameShort{chnum});
    set(gca,'FontSize',8)
    
   colormap parula

        % subplot(p(1),p(2),count+2*p(2))
    subplot(p(1),p(2),count+p(2))

hold on
    % rrn=size(dd,1);
    rrn=5;
    col=[(chnum/6)*ones(rrn,1) (1-chnum/12)*ones(rrn,1) (1-chnum/12)*ones(rrn,1)];
    rr=randperm(size(dd,1),rrn);
    for i=1:rrn
    plot(ftime(1,:),dd(rr(i),:),'LineWidth',1.5)%,"Color",col(i,:)
    end

    % p50 = mean(dd,'omitnan');
    % p25 = prctile(dd,25);
    % p75 = prctile(dd,75);
    % % % 
    % hIQR = fill([ftime(1,:),fliplr(ftime(1,:))],[p75,fliplr(p25)],"r","FaceAlpha",.1);
    % plot(ftime(1,:),p50,'LineWidth',3,"Color","r")
    set(gca, 'Ylim',clims);
    % set(gca, 'Ylim',[-.1 8]);
    set(gca, 'Ylim',[-1 13]);
    
    % pbaspect([1 1 1])
    ylabel("NFkB Dynamics")
    hold off
    set(gca, 'Xlim', [tl(1) tl(2)]);
    xticks(0:120:tl(2))
    xlabel('Time (min)')
    % ylabel(tit)
    title(gnameShort{chnum});
    set(gca,'FontSize',8)
    xlim([tl(1) tl(2)])

end

end