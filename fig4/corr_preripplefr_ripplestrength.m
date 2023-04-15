clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

basewin = [-3 -2];
win = [-1 0];
rippletype = {'medial';'lateral';'global'};

binsize = 0.01;

[r_strength,p_strength,r_length,p_length] = deal(nan(nS,3,5));
[r_strength_rpband,r_length_rpband] = deal(nan(nS,3));

for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_ripple_modulated.mat'],'fr_ripple');
    load([sdir(sessionList(iS)),'_cellTable.mat'],'T');
    load([sdir(sessionList(iS)),'_ripples.mat'],'CA1_ripple_classified','filtered_lfp');
    
    targets = tag.info.unit_id((tag.area.vis&tag.celltype.rs) | tag.area.thalamus);
    in = ismember(T.unit_id,targets);
    
    unitid = T.unit_id(in);
    spiketime = T.spike_time(in);
    
    cidx = zeros(sum(in),1);
    in = ismember(unitid,tag.info.unit_id(tag.area.thalamus));
    cidx(in) = 4; % thalamus
    [in,idx] = ismember(unitid,unit_id.vis);
    cidx(in) = cluster_idx.vis{2}(idx(in));
    
    clstnum = cellfun(@(x) sum(cidx==x),num2cell(0:4));
    if min(clstnum)<5
        continue;
    end
    
    [nss,nsstime] = deal(cell(2,1));
    for iProbe = 1:2
        if iscell(filtered_lfp.time{iProbe})
            nss_tmp = cell(length(filtered_lfp.time{iProbe}),1);
            for iw = 1:length(filtered_lfp.time{iProbe})
                nss_tmp{iw} =  mean(NormalizedSquaredSignal_HJ([filtered_lfp.time{iProbe}{iw}',...
                    filtered_lfp.lfp{iProbe}{iw}]),2);
            end
            nss{iProbe} = cell2mat(nss_tmp);
            nsstime{iProbe} = cell2mat(filtered_lfp.time{iProbe}');
        else
            nss{iProbe} = mean(NormalizedSquaredSignal_HJ([filtered_lfp.time{iProbe}',...
                filtered_lfp.lfp{iProbe}]),2);
            nsstime{iProbe} = filtered_lfp.time{iProbe};
        end
    end
        
    for iRp = 1:3
        preripplespk = cellfun(@(y) cellfun(@(x) sum(y>=x+win(1) & y<=x+win(2)),...
            num2cell(CA1_ripple_classified.(rippletype{iRp})(:,1)))',...
            spiketime,'UniformOutput',false);
        baseripplespk = cellfun(@(y) cellfun(@(x) sum(y>=x+basewin(1) & y<=x+basewin(2)),...
            num2cell(CA1_ripple_classified.(rippletype{iRp})(:,1)))',...
            spiketime,'UniformOutput',false);
        preripplespkz = cell2mat(cellfun(@(x,y) (x-mean(y))/std(y),...
            preripplespk,baseripplespk,'UniformOutput',false));
        preripplespkz(isinf(preripplespkz)) = nan;

    
        if iRp==3
            ripplestrength = max(CA1_ripple_classified.global_strength,[],2);
            ripplelength = diff(CA1_ripple_classified.global,[],2);
            
            [preripplenss,baseripplenss] = deal(nan(size(CA1_ripple_classified.(rippletype{iRp}),1),2));
            for iProbe = 1:2
                intime = cellfun(@(x) nsstime{iProbe}>=x+win(1)-0.1 &...
                    nsstime{iProbe}<=x+win(2)+0.1,...
                    num2cell(CA1_ripple_classified.(rippletype{iRp})(:,1)),'UniformOutput',false);
                preripplenss(:,iProbe) = cell2mat(cellfun(@(x,y) sum(interp1(nsstime{iProbe}(y),...
                    nss{iProbe}(y),x+[win(1)+binsize/2:binsize:win(2)-binsize/2]))/(1/binsize),...
                    num2cell(CA1_ripple_classified.(rippletype{iRp})(:,1)),intime,'UniformOutput',false));
                inbasetime = cellfun(@(x) nsstime{iProbe}>=x+basewin(1)-0.1 &...
                    nsstime{iProbe}<=x+basewin(2)+0.1,...
                    num2cell(CA1_ripple_classified.(rippletype{iRp})(:,1)),'UniformOutput',false);
                baseripplenss(cellfun(@sum,inbasetime)~=0,iProbe) = cell2mat(cellfun(@(x,y) sum(interp1(nsstime{iProbe}(y),...
                    nss{iProbe}(y),x+[basewin(1)+binsize/2:binsize:basewin(2)-binsize/2]))/(1/binsize),...
                    num2cell(CA1_ripple_classified.(rippletype{iRp})(cellfun(@sum,inbasetime)~=0,1)),inbasetime(cellfun(@sum,inbasetime)~=0),'UniformOutput',false));
            end
            preripplenss = max(preripplenss,[],2);
%             preripplenss = (mean(preripplenss)-nanmean(mean(baseripplenss,2)))/nanstd(mean(baseripplenss,2));
        else
            ripplestrength = CA1_ripple_classified.(rippletype{iRp})(:,4);
            ripplelength = diff(CA1_ripple_classified.(rippletype{iRp})(:,[1,3]),[],2);
            
            intime = cellfun(@(x) nsstime{iRp}>=x+win(1)-0.1 &...
                nsstime{iRp}<=x+win(2)+0.1,...
                num2cell(CA1_ripple_classified.(rippletype{iRp})(:,1)),'UniformOutput',false);
            preripplenss = cell2mat(cellfun(@(x,y) sum(interp1(nsstime{iRp}(y),...
                nss{iRp}(y),x+[win(1)+binsize/2:binsize:win(2)-binsize/2]))/(1/binsize),...
                num2cell(CA1_ripple_classified.(rippletype{iRp})(:,1)),intime,'UniformOutput',false));
            inbasetime = cellfun(@(x) nsstime{iRp}>=x+basewin(1)-0.1 &...
                nsstime{iRp}<=x+basewin(2)+0.1,...
                num2cell(CA1_ripple_classified.(rippletype{iRp})(:,1)),'UniformOutput',false);
            baseripplenss = nan(size(CA1_ripple_classified.(rippletype{iRp}),1),1);
            baseripplenss(cellfun(@sum,inbasetime)~=0) = cell2mat(cellfun(@(x,y) sum(interp1(nsstime{iRp}(y),...
                nss{iRp}(y),x+[basewin(1)+binsize/2:binsize:basewin(2)-binsize/2]))/(1/binsize),...
                num2cell(CA1_ripple_classified.(rippletype{iRp})(cellfun(@sum,inbasetime)~=0,1)),...
                inbasetime(cellfun(@sum,inbasetime)~=0,1),'UniformOutput',false));
%             preripplenss = (preripplenss-nanmean(baseripplenss))/nanstd(baseripplenss);
        end
        
        for iclst = 1:5
            [r_strength(iS,iRp,iclst),p_strength(iS,iRp,iclst)] =...
                corr(nanmean(preripplespkz(cidx==iclst-1,:))',ripplestrength,'rows','complete');
            [r_length(iS,iRp,iclst),p_length(iS,iRp,iclst)] =...
                corr(nanmean(preripplespkz(cidx==iclst-1,:))',ripplelength,'rows','complete');
        end
        
        r_strength_rpband(iS,iRp) = corr(preripplenss,ripplestrength,'rows','complete');
        r_length_rpband(iS,iRp) = corr(preripplenss,ripplelength,'rows','complete');
    end
end

%%
dir = 'D:\OneDrive - University of California, San Francisco\figures\2.allen';
pttest = nan(3,4);
tttest = nan(3,4);
for iRp = 1:3
    if iRp<3
        if iRp==1
            fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5 3.5]);
        end
        axes('Position',axpt(2,1,iRp,1,axpt(10,10,2:10,1:9)));
    elseif iRp==3
        fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3.75 3.5]);
    end
    rdata = squeeze(r_strength(:,iRp,[4,5,2,3]));
    hold on;
    plot(1:4,rdata,'color',[0.6 0.6 0.6],'LineWidth',0.35);
    errorbar(1:4,nanmean(rdata),nanstd(rdata)/sqrt(sum(~isnan(rdata(:,1)))),'k','CapSize',3);
    [h,pttest(iRp,:),~,stat] = ttest(rdata);
    tttest(iRp,:) = stat.tstat;
    for iclst = 1:4
       if h(iclst)==1
           text(iclst,0.5,'*','FontSize',8);
       end
    end
    plot([0 6],[0 0],'k:');
    ylim([-0.4 0.4])
    xlim([0.5 4.5])
    set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',...
        1:4,'XTickLabel',{'Inh','Thal','iAct','dAct'},'XTickLabelRotation',45,'YTick',-0.4:0.4:0.4);
    if rem(iRp,2)==1
        ylabel({'Correlation'; '(pre-ripple firing'; 'vs. ripple strength)'});
    else
        set(gca,'YTickLabel',[]);
        print(fHandle,'-depsc',[dir,'\corr_prespk_ripplestrength_localizedripples.ai']);
    end
end
print(fHandle,'-depsc',[dir,'\fig4\corr_prespk_ripplestrength_globalripples.ai']);
% 
%%
iw = 1;
for iRp = 1:3
    if iRp<3
        if iRp==1
            fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4.5 3.5]);
        end
        subplot(1,2,iRp)
    elseif iRp==3
        fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 3.5]);
    end
    rdata = squeeze(r_strength_rpband(:,iRp));
    hold on;
    bar(0.5,nanmean(rdata),'FaceColor',[0.6 0.6 0.6],'EdgeColor','k');
    scatter(rand(size(rdata,1),1)*0.6+0.2,rdata,2,'k','filled');
    errorbar(0.5,nanmean(rdata),nanstd(rdata)/sqrt(sum(~isnan(rdata(:,1)))),'k','CapSize',3);
    h = ttest(rdata);
    if h==1
        text(0.5,0.5,'*','FontSize',8);
    end
%     plot([0 6],[0 0],'k:');
    ylim([-0.4 0.4])
    xlim([-0.5 1.5])
    set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
        'YTick',-0.4:0.4:0.4,'XTick',[]);
    if rem(iRp,2)==1
        ylabel({'Correlation'; '(pre-ripple ripple power'; 'vs. ripple strength)'});
    else
        set(gca,'YTickLabel',[]);
        print(fHandle,'-depsc',[dir,'\corr_preripple_ripplestrength_localizedripples_',num2str(iw),'.ai']);
    end
end
print(fHandle,'-depsc',[dir,'\corr_preripple_ripplestrength_globalripples_',num2str(iw),'.ai']);

