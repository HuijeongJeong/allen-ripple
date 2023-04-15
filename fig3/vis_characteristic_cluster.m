clc; clearvars; close all;

rng(3);

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
Tfst = readtable('D:\OneDrive\1.allen-andermann\time_to_first_spike.csv');

refList = {'medial';'lateral'};
typeList = {'rs';'fs'};

nclst = 3;

[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
session_id = tag.info.session_id(idx);

unitid = unit_id.vis;
sessionList = unique(session_id);
nS = length(sessionList);

%%
[pval_rf,ctype,cluster,area,tfs,rfsize,modidx,fr_dg] = deal(cell(nS,1));

for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    
    invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis));
    inthal = ismember(T.unit_id,tag.info.unit_id(ismember(tag.info.structure,{'LGd';'LGv';'LP'})));
    
    [in,idx] = ismember(T.unit_id,unitid);
    clstidx = NaN(length(T.unit_id),1);
    clstidx(in) = cluster_idx.vis{nclst-1}(idx(in));
    
    [~,idx] = ismember(T.unit_id,tag.info.unit_id);
    celltype = [tag.celltype.rs(idx), tag.celltype.fs(idx)];
    
    in = invis;
  
    ctype{iS} = celltype(in,:);
    cluster{iS} = clstidx(in);
    area{iS} = T.ecephys_structure_acronym(in);
    pval_rf{iS} = T.p_value_rf(in);
    fr_dg{iS} = T.firing_rate_dg(in);
        
    [~,idx] = ismember(T.unit_id(in),Tfst.ecephys_unit_id);
    tfs{iS} = Tfst.time_to_first_spike_fl(idx);
    rfsize{iS} = T.area_rf(in);
    modidx{iS} = T.mod_idx_dg(in);
end

%%

ctype = cell2mat(ctype);
cluster = cell2mat(cluster);
cluster(isnan(cluster)) = 0;

area = cat(1,area{:});
pval_rf = cell2mat(pval_rf);
fr_dg = cell2mat(fr_dg);
tfs = cell2mat(tfs);
rfsize = cell2mat(rfsize);
modidx = cell2mat(modidx);

nCell = 100;
nIter = 100;

%%
close all
measureList = {'Receptive field area (deg^2)','Time to first spike (ms)','log_1_0 MI'};
ct = [0.6 0.6 0.6; cbrewer('qual','Dark2',3)];
data = [rfsize,tfs*1000,log10(modidx)];
binrange = {0:100:2500; 30:5:100; -3:0.1:1.5};
iCT = 1;

sigfon = NaN(nIter,nclst+1);
data_cumsum = cell(3,nclst+1);

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20/2 20/6]);
h(1) = axes('Position',axpt(3,1,1,1,axpt(1,10,1,1:9))); hold on;
h(2) = axes('Position',axpt(3,1,2,1,axpt(1,10,1,1:9))); hold on;
h(3) = axes('Position',axpt(3,1,3,1,axpt(1,10,1,1:9))); hold on;
for iClst = 1:nclst+1
    in = find(cluster==iClst-1 & ctype(:,iCT));
    
    for iI = 1:nIter
        initer = randsample(in,nCell);
        pvaliter = pval_rf(initer);
        dataiter = data(initer,:);
        
        sigvis = pvaliter<0.01 & dataiter(:,1)<2500 & dataiter(:,2)<=100;
        sigfon(iI,iClst) = sum(sigvis)/nCell*100;
        for id = 1:3
            if isempty(data_cumsum{id,iClst})
                data_cumsum{id,iClst} = NaN(nIter,length(binrange{id}));
            end
            data_cumsum{id,iClst}(iI,:) = cumsum(histc(dataiter(sigvis,id),binrange{id}))/sum(sigvis)*100;
        end
    end
    
    for id = 1:3
        m = mean(data_cumsum{id,iClst});
        data_cumsum_sort = sort(data_cumsum{id,iClst});
        s = data_cumsum_sort([0.05 0.95]*nIter,:);
        fill([binrange{id}, flip(binrange{id})],[s(1,:) flip(s(2,:))],ct(iClst,:),...
            'EdgeColor','none','Parent',h(id));
        plot(binrange{id},m,'Color',ct(iClst,:),'Parent',h(id));
        if iClst==nclst+1
           alpha(h(id),0.2);
           set(h(id),'XLim',[min(binrange{id}) max(binrange{id})]);
           xlabel(h(id),measureList{id});
        end
    end
end
set(h,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'YTick',0:20:100);
set(h(2:3),'YTickLabel',[])
ylabel(h(1),'Cumulative fraction (%)');
cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig3');
print(fHandle,'-depsc','-painters','visual_characteristics.ai');

group = repmat([1 2 3 4],nIter,1);
[p,~,stats] = anova1(sigfon(:),group(:));
pmulti = multcompare(stats);

%%
n = nan(4,2);
for iClst = 1:4
inclst = cluster==iClst-1 & ctype(:,iCT);
sig = pval_rf<0.01 & rfsize<2500 & tfs*1000<=100;
n(iClst,1) = sum(inclst & ~sig);
n(iClst,2) = sum(inclst & sig);
end
[p,chi2stat,proportionT] = chisquare(n);

p_sub = nan(4,4);
for iclst = 1:3
    for jclst = iclst+1:4 
p_sub(iclst,jclst) = chisquare(n([iclst,jclst],:)')*6;
    end
end

%%
sig = pval_rf<0.01 & rfsize<2500 & tfs*1000<=100;
rftest = rfsize(sig & ctype(:,iCT));
clustertest = cluster(sig & ctype(:,iCT));
anova1(rftest,clustertest);

%%
close all
m = mean(sigfon);
sigfon_sort = sort(sigfon);
s = sigfon_sort([0.05 0.95]*nIter,:);

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 20*0.19]); 
hold on;
for iClst = 1:nclst+1
   plot([iClst iClst],s(:,iClst),'Color',ct(iClst,:));
   scatter(iClst,m(iClst),2,ct(iClst,:),'s','filled');
end
ylim([0 60])
xlim([0.5 4.5]);
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',1:4,...
    'XTickLabel',{'Nomod','Cluster 1','Cluster 2','Cluster 3'},'XTickLabelRotation',45,...
    'YTick',0:20:60);
ylabel('FON (%)');
cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig3');
print(fHandle,'-depsc','-painters','sig_visual_nurons.ai');
%%
close all
nCell = 30;
areaList = {'VISp';'VISl';'VISrl';'VISal';'VISpm';'VISam'};
nArea = length(areaList);
ylimit = {[200 1200];[48 92]; [-0.6 0.9]};
yticks = {200:200:1200;50:10:90;-0.6:0.3:0.9};

data_median = cell(3,nclst+1);
x = [repmat(1:nArea,2,1);NaN(1,nArea)];
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.5 20/6]);
h(1) = axes('Position',axpt(3,1,1,1,axpt(1,10,1,1:9),[0.13 0.05])); hold on;
h(2) = axes('Position',axpt(3,1,2,1,axpt(1,10,1,1:9),[0.13 0.05])); hold on;
h(3) = axes('Position',axpt(3,1,3,1,axpt(1,10,1,1:9),[0.13 0.05])); hold on;
for iClst = 1:nclst+1
    for iA = 1:nArea
        in = find(cluster==iClst-1 & ctype(:,iCT) & strcmp(area,areaList{iA}));
        for iI = 1:nIter
            initer = randsample(in,nCell);
            pvaliter = pval_rf(initer);
            dataiter = data(initer,:);
            
            sigvis = pvaliter<0.01 & dataiter(:,1)<2500 & dataiter(:,2)<=100;
            sigfon(iI,iA) = sum(sigvis)/nCell*100;
            for id = 1:3
                if isempty(data_median{id,iClst})
                   data_median{id,iClst} = NaN(nIter,nArea);
                end
                data_median{id,iClst}(iI,iA) = median(dataiter(sigvis,id));
            end
        end
    end
    
    m = cellfun(@nanmean,data_median(:,iClst),'UniformOutput',false);
    s = cellfun(@sort,data_median(:,iClst),'UniformOutput',false);
    s = cellfun(@(x) [x([0.05 0.95]*nIter,:);NaN(1,nArea)],s,'UniformOutput',false);
     
    for id = 1:3
        plot([1:nArea]+(iClst-1)*0.15,m{id},'Color',ct(iClst,:),'Parent',h(id));
        plot([x(:)]+(iClst-1)*0.15,s{id}(:),'Color',ct(iClst,:),'Parent',h(id));
        if iClst==nclst+1
           xlim(h(id),[0.3 nArea+1]);
           set(h(id),'YLim',ylimit{id},'YTick',yticks{id});
           ylabel(h(id),measureList{id});
        end
    end
end
set(h,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'Box','off',...
    'XTick',[1:6]+0.18,'XTickLabel',{'V1';'LM';'RL';'AL';'PM';'AM'},'XTickLabelRotation',45)
cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig3');
% print(fHandle,'-depsc','-painters','visual_characteristics_area.ai');

%%
data = cellfun(@(x) x(:,[1,6]),data_median(id,:),'UniformOutput',false);
simple_mixed_anova(cat(3,data_median{id,:}))








