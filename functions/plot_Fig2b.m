function plot_Fig2b(data,gname,tl,clims,p,cc,bb)
rng(6)
cString= {
    'Sum of P65 in Nuc',...
    'Mean of P65 in Nuc',...
    'Median of P65 in Nuc',...
    'Sum of P65 in Around Nuc',...
    'Mean of P65 in Around Nuc',...
    'Median of P65 in Around Nuc',...
    'Sum of P65 in Cytosol',...
    'Mean of P65 in Cytosol',...
    'Median of P65 in Cytosol',...
    'Mean of P65 around Center',...
    'Median of P65 around Center',...
    'Std of P65 around Center ',...
    'Total P65',...
    'Sum Nuc/Sum Cyt',...
    'Mean/Mean',...
    'Median/Median',...
    'Mean/Mean Around Nuc',...
    'Median/Median Around Nuc'
    };

tit=cString{cc};

figure("Name",tit)
count=0;


for chnum=1:numel(gname)
    if mod(count+1,17)==0
        if b==0
        h = axes('visible','off'); 
        c = colorbar(h,'Position',[0.93 0.12 0.022 0.7]);  % attach colorbar to h
        caxis(h,clims); 
        end
        
        figure
        count=0;
    end
count=count+1;

dd=[data.("NfkB") ];
ftime=[data.("Time")]; 
zz=data.Category==gname{chnum};
ftime=ftime(zz,:);
dd=dd(zz,:);
ncells=size(dd,1);


% auc=trapz(dd,2);
auc=trapz(dd(:,1:20),2);
[~,Ind]=sort(auc,'descend');

if bb==0
    subplot(p(1),p(2),count)
    imagesc(ftime(1,:),1:ncells,dd(Ind,:),clims)
    set(gca, 'Ylim', [0 ncells]);
    yticks([])
    ylabel([])
    % pbaspect([1 1 1])
    set(gca, 'Xlim', [tl(1) tl(2)]);
    xticks(0:120:tl(2))
    xlabel('Time (min)')
    % % ylabel(tit)
    % title(gname{chnum});
    % set(gca,'FontSize',10,'FontName','Times New Roman')
    colormap parula

else
    subplot(p(1),p(2),count)
    hold on
    % rrn=size(dd,1);
    % rrn=5;
    % col=[(ceil(chnum/4)/4)*ones(rrn,1) (1-ceil(chnum/4)/4)*ones(rrn,1) (1-ceil(chnum/4)/4)*ones(rrn,1)];
    % 
    % rr=randperm(size(dd,1),rrn);
    % % plot(ftime(1,:),dd(rr,:)','LineWidth',1,'Color',col)
    % for i=1:rrn
    % plot(ftime(1,:),dd(rr(i),:),'LineWidth',1.5)%,"Color",col(i,:)
    % end

    p50 = median(dd,'omitnan');
    % p50=p50-min(p50);
    % stdNf = .5*std(dd,'omitnan');
    % p25 = p50-stdNf;
    % p75 = p50+stdNf;
    p25 = prctile(dd,25);
    p75 = prctile(dd,75);

    hIQR = fill([ftime(1,:),fliplr(ftime(1,:))],[p75,fliplr(p25)],"r","FaceAlpha",.1);
    plot(ftime(1,:),p50,'LineWidth',3,'Color','r')

    if chnum<5
        % ylim([-1 6])
        % yticks([0 3 6])
        ylim([-.1 4])
    elseif chnum>=5 && chnum<9
        % ylim([-1 12])
        % yticks([0 6 12])
        ylim([-.1 8])
    else
        % ylim([-1 10])
        % yticks([0 5 10])
        ylim([-.1 6])
    end

    % pbaspect([1 1 1])
    % ylabel("NFkB Dynamics")
    hold off
    set(gca, 'Xlim', [tl(1) tl(2)]);
    xticks(0:120:tl(2))
    xlabel('Time (min)')
    % ylabel(tit)
    % title(gname{chnum});
    set(gca,'FontSize',10,'FontName','Times New Roman')

end

end

end