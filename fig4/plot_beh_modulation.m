clc; clearvars; close all;
rng(3);

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

refList = {'medial';'lateral'};
typeList = {'rs';'fs'};

nclst = 3;

[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
session_id = tag.info.session_id(idx);

unitid = unit_id.vis;
sessionList = unique(session_id);
nS = length(sessionList);

%%
[behavior_corr,behavior_corrp,ctype,cluster,area,...
    behavior_selec,behavior_selecp,laminar_loc] = deal(cell(nS,1));
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_stimulus_evoked.mat'],'fr_behavior','fr_movie');
    load([sdir(sessionList(iS)),'_Events.mat'],'drifting_gratings','stimulus_presentation');
    load([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_win');
    
    invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis));
    
    [in,idx] = ismember(T.unit_id,unitid);
    clstidx = NaN(length(T.unit_id),1);
    clstidx(in) = cluster_idx.vis{nclst-1}(idx(in));
    
    [~,idx] = ismember(T.unit_id,tag.info.unit_id);
    celltype = [tag.celltype.rs(idx), tag.celltype.fs(idx)];
    laminar_loc{iS} = tag.info.distance_from_L4(idx);

    laminar_loc{iS} = laminar_loc{iS}(invis);
    ctype{iS} = celltype(invis,:);
    cluster{iS} = clstidx(invis);
    area{iS} = T.ecephys_structure_acronym(invis);
  
    behavior_corr{iS} = [fr_behavior.run.corr.r(invis), fr_behavior.pupil.corr.r(invis),...
        fr_behavior.pupil_norun.corr.r(invis)];
    behavior_corrp{iS} = [fr_behavior.run.corr.pval(invis), fr_behavior.pupil.corr.pval(invis),...
        fr_behavior.pupil_norun.corr.pval(invis)];
    
    behavior_selec{iS} = [fr_behavior.run.selec.idx(invis), fr_behavior.pupil.selec.idx(invis),...
        fr_behavior.pupil_norun.selec.idx(invis)];
    behavior_selecp{iS} = [fr_behavior.run.selec.pval(invis), fr_behavior.pupil.selec.pval(invis),...
        fr_behavior.pupil_norun.selec.pval(invis)];
end

%%
cluster = cell2mat(cluster);
cluster(isnan(cluster)) = 0;
celltype = cell2mat(ctype);
area = cat(1,area{:});
behavior_selec = cell2mat(behavior_selec);
behavior_selecp = cell2mat(behavior_selecp);

%%
ct = [0.6 0.6 0.6; cbrewer('qual','Dark2',3)];
iCT = 1;
nCell = 100;
nIter = 100;
binrange = -1.5:0.05:1.5;
titleList = {'Running';'Pupil'};
%%

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20/2 20*0.14]);
h(1) = axes('Position',axpt(2,1,1,1,axpt(3,10,1:2,1:9,[],[0.1 0.05]),[0.07 0.05])); hold on;
h(2) = axes('Position',axpt(2,1,2,1,axpt(3,10,1:2,1:9,[],[0.1 0.05]),[0.07 0.05])); hold on;
h(3) = axes('Position',axpt(2,1,1,1,axpt(3,10,3,1:9,[],[0.1 0.05]),[0.1 0.05])); hold on;
h(4) = axes('Position',axpt(2,1,2,1,axpt(3,10,3,1:9,[],[0.1 0.05]),[0.1 0.05])); hold on;
for iClst = 1:4
    in = find(cluster==iClst-1 & celltype(:,iCT));
    beh_cumsum = NaN(nIter,length(binrange),2);
    beh_sig = NaN(nIter,2);
    for iIter = 1:nIter
        initer = randsample(in,nCell);
        selec_sub = behavior_selec(initer,:);
        selecp_sub = behavior_selecp(initer,:);
        beh_cumsum(iIter,:,1) = cumsum(histc(selec_sub(:,1),binrange))/nCell; % running
        beh_sig(iIter,1) = mean(selecp_sub(:,1)<0.05);
        beh_cumsum(iIter,:,2) = cumsum(histc(selec_sub(:,3),binrange))/nCell; % pupil - no running
        beh_sig(iIter,2) = mean(selecp_sub(:,3)<0.05);
    end
 
    for i = 1:2
        m = mean(squeeze(beh_cumsum(:,:,i)*100));
        beh_cumsum_sort = sort(squeeze(beh_cumsum(:,:,i)*100));
        s = beh_cumsum_sort([0.05 0.95]*nIter,:);
        fill([binrange flip(binrange)],[s(1,:) flip(s(2,:))],ct(iClst,:),...
            'EdgeColor','none','Parent',h(i));
        plot(binrange,m,'Color',ct(iClst,:),'Parent',h(i));
        
        m = mean(beh_sig(:,i)*100);
        behsig_sort = sort(beh_sig(:,i)*100);
        s = behsig_sort([0.05 0.95]*nIter);
        plot([iClst iClst],s,'Color',ct(iClst,:),'Parent',h(i+2));
        scatter(iClst,m,2,ct(iClst,:),'s','filled','Parent',h(i+2));
        if iClst==4
            plot([0 0],[0 100],'k:','LineWidth',0.35,'Parent',h(i));
            title(h([i,i+2]),titleList{i});
        end
    end
end
set(h,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'YTick',0:20:100);
set(h(1),'XTick',-1:0.5:1,'XLim',[-1.2 1.2]);
set(h(2),'XTick',-0.6:0.3:0.6,'XLim',[-0.6 0.6]);
set(h(3:4),'XTick',1:4,'XTickLabel',{'Nomod';'Cluster 1';'Cluster 2';'Cluster 3'},...
    'XTickLabelRotation',45,'YTick',0:20:100,'YLim',[0 100],'XLim',[0.5 4.5]);
set(h([2,4]),'YTickLabel',[]);
xlabel(h(1:2),'Selectivity index');
ylabel(h(1),'Cumulative fraction (%)');
ylabel(h(3),'FON (%)');
alpha(h(1),0.2);
alpha(h(2),0.2);
% cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig2')
% print(fHandle,'-depsc','-painters','behavior_cumfraction.ai');

%%
close all
nCell = 30;
areaList = {'VISp';'VISl';'VISrl';'VISal';'VISpm';'VISam'};
nArea = length(areaList);
xx = [repmat(1:nArea,2,1);NaN(1,nArea)];
beh_sig = cell(nArea,1);

% beh_median = NaN(nIter,nArea,2);
clear h
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.26 20*0.16]);
h(1) = axes('Position',axpt(2,1,1,1,[],[0.1 0.05])); hold on;
h(2) = axes('Position',axpt(2,1,2,1,[],[0.1 0.05])); hold on;
for iClst = 1:4
    for iA = 1:nArea
        in = find(cluster==iClst-1 & celltype(:,iCT) & strcmp(area,areaList{iA}));
        
        beh_sig{iA} = behavior_selec(in,:);
%         for iIter = 1:nIter
%             initer = randsample(in,nCell);
%             selec_sub = behavior_selec(initer,:);
%             beh_median(iIter,iA,1) = nanmedian(selec_sub(:,1)); % running
%             beh_median(iIter,iA,2) = nanmedian(selec_sub(:,3)); % pupil - no running
%         end
    end
    for i = [1,3]
        m = cellfun(@(x) nanmean(x(:,i)),beh_sig);
        s = cellfun(@(x) nanstd(x(:,i))/sqrt(size(x,1)),beh_sig);
%         m = mean(squeeze(beh_sig(:,:,i)));
%         beh_median_sort = sort(squeeze(beh_median(:,:,i)));
%         s = [beh_median_sort([0.05 0.95]*nIter,:);NaN(1,nArea)];
        
%         m = mean(squeeze(beh_median(:,:,i)));
%         beh_median_sort = sort(squeeze(beh_median(:,:,i)));
%         s = [beh_median_sort([0.05 0.95]*nIter,:);NaN(1,nArea)];
        errorbar([1:6]+(iClst-1)*0.12,m,s,'Color',ct(iClst,:),'CapSize',0,'Parent',h(round(i/2)));
%         plot(x(:),s(:),'Color',ct(iClst,:),'Parent',h(i));
%         plot([1:nArea]+(iClst-1)*0.12,m,'Color',ct(iClst,:),'Parent',h(i));
        %         scatter([1:nArea]+(iClst-1)*0.12,m,3,ct(iClst,:),'s','filled','Parent',h(i));
        if iClst==4
            title(h(round(i/2)),titleList{round(i/2)});
            plot([0.3 nArea+1],[0 0],'k:','Parent',h(round(i/2)));
        end
    end
    
end
set(h,'XLim',[0.36 nArea+1],'XTick',[1:6]+0.18,'XTickLabel',{'V1';'LM';'RL';'AL';'PM';'AM'},...
    'XTickLabelRotation',45,'Box','off','TickDir','out','FontSize',7);
set(h(1),'YTick',-0.6:0.3:0.6,'YLim',[-0.6 0.6]);
set(h(2),'YTick',-0.2:0.1:0.2,'YLim',[-0.2 0.2]);
ylabel(h(1),'Selectivity index');
% cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig2')
% print(fHandle,'-depsc','-painters','behavior_area.ai');
%%
in = celltype(:,iCT) & ismember(area,areaList);
[~,table,stat] = anovan(behavior_selec(in,3),{cluster(in),area(in)},'model','interaction');
multcompare(stat)
%%
binrange = -1:0.05:1;
ct = [0 0 0; cbrewer('qual','Dark2',3)];
areaList = {'VISp';'VISl';'VISrl';'VISal';'VISpm';'VISam'};
titleList = {'running';'pupil';'pupil w/o running'};
iCT = 1;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 12 7]);
for iBeh = 1:3
    h(iBeh) = axes('Position',axpt(3,2,iBeh,1,[],[0.1 0.15])); hold on;
    h(iBeh+3) = axes('Position',axpt(3,2,iBeh,2,[],[0.1 0.2])); hold on;
end
for iClst = 1:4
    if iClst==1
        in = cellfun(@(x,y) isnan(x) & y(:,iCT),...
            cluster,ctype,'UniformOutput',false);
    else
        in = cellfun(@(x,y) x==iClst-1 & y(:,iCT),...
            cluster,ctype,'UniformOutput',false);
    end
    data_p = cell2mat(cellfun(@(x,y) x(y,:),behavior_selecp,in,'UniformOutput',false));
    data_r = cell2mat(cellfun(@(x,y) x(y,:),behavior_selec,in,'UniformOutput',false));
    
    [m,s] = deal(NaN(length(areaList),3));
    for iArea = 1:length(areaList)
        inarea = cellfun(@(x,y) strcmp(x,areaList{iArea}) & y,...
            area,in,'UniformOutput',false);
        dataarea = cell2mat(cellfun(@(x,y) x(y,:),behavior_selec,inarea,...
            'UniformOutput',false));
        m(iArea,:) = nanmean(dataarea);
        s(iArea,:) = nanstd(dataarea)/sqrt(length(dataarea));
    end
    for iBeh = 1:3
       bincount = histc(data_r(:,iBeh),binrange);
       plot(binrange,cumsum(bincount)/sum(bincount),'Color',ct(iClst,:),'Parent',h(iBeh));
       errorbar(1:length(areaList),m(:,iBeh),s(:,iBeh),'CapSize',3,'Color',ct(iClst,:),'Parent',h(iBeh+3));
       if iClst==4
           plot([0 0],[0 1],'k:','Parent',h(iBeh));
           title(titleList{iBeh},'Parent',h(iBeh));
       end
       xlabel(h(iBeh),'Modulation index');
    end
end
ylabel(h(4),'Modulation index');
ylabel(h(1),'Cumulative fraction');
set(h,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,'XTick',-1:1:1);
set(h(1:3),'XLim',[-1 1]);
set(h(4:6),'XLim',[0 length(areaList)+1],'XTick',1:6,'XTickLabel',areaList,'XTickLabelRotation',45);

% cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig2');
% print(fHandle,'-dtiff','-r600','behavior_modulation.tif');
