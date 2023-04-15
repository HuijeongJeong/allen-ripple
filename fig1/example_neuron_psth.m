clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');

%% setting up 
k = 3; % number of cluster 

[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
celltype = [tag.celltype.rs(idx),tag.celltype.fs(idx)];
sessionid = tag.info.session_id(idx);
sessionList = unique(sessionid);
nS = length(sessionList);

win = [-2 2];
bin = 0.01;
resolution = 10;
ct = {'r';'b'};
%%
iCT = 1;
for iS = 5:9
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_ripples.mat'],'CA1_ripple_classified');
    
    [in,idx] = ismember(T.unit_id,unit_id.vis(celltype(:,iCT)));
    cidx = cluster_idx.vis{k-1}(idx(in));
    spiketime = T.spike_time(in);
    unitid = T.unit_id(in);

    [cidx,sortidx] = sort(cidx);
    spiketime = spiketime(sortidx);
    unitid = unitid(sortidx);

    rippletime = [CA1_ripple_classified.medial(:,1);CA1_ripple_classified.lateral(:,1);...
        CA1_ripple_classified.global(:,1)];
    rippleidx = [ones(size(CA1_ripple_classified.medial,1),1);...
        ones(size(CA1_ripple_classified.lateral,1),1)*2;...
        ones(size(CA1_ripple_classified.global,1),1)*3];
    rippleidx = [rippleidx==1, rippleidx==2, rippleidx==3];
    

    fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20 20]);
    for iunit = 1:sum(in)
        spikeripple = spikeWin(spiketime{iunit},rippletime,win);
        [xpt,ypt,time,~,spkconv,~,spksem] = rasterPSTH(spikeripple,rippleidx,win,bin,resolution);

        x = rem(iunit,10);
        if x==0
            x = 10;
        end
        h(1) = axes('Position',axpt(1,2,1,1,axpt(10,10,x,ceil(iunit/10),[],[0.025 0.025])));hold on;
        h(2) = axes('Position',axpt(1,2,1,2,axpt(10,10,x,ceil(iunit/10),[],[0.025 0.025])));hold on;
        for iref = 1:2
            scatter(xpt{iref},ypt{iref},0.5,'.',ct{iref},'Parent',h(1));
            fill([time flip(time)],spksem(iref,:),ct{iref},'EdgeColor','none','Parent',h(2));
            plot(time,spkconv(iref,:),'Color',ct{iref},'Parent',h(2));
        end
        y2 = get(h(2),'YLim');
        if y2(2)<2
            y2(2) = 2;
        end
        plot([0 0],y2,'k:','Parent',h(2));
        plot([0 0],[1 find(rippleidx(:,2),1,'last')],'k:','Parent',h(1))
        set(h,'Box','off','TickDir','out','XLim',[-1.5 1.5],'FontSize',4)
        set(h(1),'YTick',find(rippleidx(:,2),1,'last'),'YLim',[1 find(rippleidx(:,2),1,'last')])
        set(h(2),'YTick',[ceil(y2(1)),floor(y2(2))],'YLim',y2)
        alpha(0.2);
        if ceil(iunit/10)==ceil(sum(in)/10) 
            xlabel(h(2),'Time (s)');
        else
            set(h(2),'XTickLabel',[])
        end
        if x==1
            ylabel(h(1),'Ripple');
            ylabel(h(2),'Rate (Hz)')
        end
    end
    print(fHandle,'-dtiff','-r600',[fileparts(sdir(sessionList(iS))),'\figures\ripple_modulated_units.tif']);
end