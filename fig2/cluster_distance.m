clearvars; clc; close all;

load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');

k = 3;
sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
celltype = [tag.celltype.rs(idx), tag.celltype.fs(idx)];
ccf = tag.info.ccf(idx,:);
session = tag.info.session_id(idx);
area = tag.info.structure(idx);
probeid = tag.info.probe_id(idx);

%%
iCT = 1;
fraction_wiarea = NaN(nS,2);
[distance_wiarea,distance_bwarea,distance_sameprobe] = deal(cell(nS,2));
for iS = 1:nS
    iS
    insession = session==sessionList(iS) & celltype(:,iCT) & sum(ccf<0,2)==0;
    nunit = sum(insession);
    
    if nunit<=5
        continue;
    end
    
    ccf_session = ccf(insession,:);
    cidx_session = cluster_idx.vis{2}(insession);
    
    dist = squareform(pdist(ccf_session));
    dist(probeid(insession)==probeid(insession)') = NaN;
    dist(logical(triu(ones(size(dist,1)),0))) = NaN;
    
    x = repmat(cidx_session,1,nunit);
    y = repmat(cidx_session',nunit,1);
    wicluster = x(:)==y(:);
    
%      fraction_wiarea(iS,:) = [mean(wiarea(wicluster)),mean(wiarea(~wicluster))];
%      distance_wiarea(iS,:) = {dist(wiarea&wicluster), dist(wiarea&~wicluster)};
%      distance_bwarea(iS,:) = {dist(~wiarea&wicluster), dist(~wiarea&~wicluster)};
     distance_sameprobe(iS,:) = {dist(wicluster), dist(~wicluster)};
end

%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 20*0.18]);
hold on;
plot(1:2,fraction_wiarea*100,'Color',[0.6 0.6 0.6],'LineWidth',0.35);
errorbar(1:2,nanmean(fraction_wiarea*100),nanstd(fraction_wiarea*100)/sqrt(sum(~isnan(fraction_wiarea(:,1)))),'k');
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',1:2,...
    'XTickLabel',{'w/i cluster';'b/w cluster'},'XTickLabelRotation',45,'YTick',0:10:40)
ylim([0 45])
xlim([0.5 2.5])
ylabel({'Fraction of';'within-area pair (%)'})
cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig2');
print(fHandle,'-depsc','-painters','withinarea_pair_fraction.ai');

%%
clr = {[1 0.6 0.6],[1 0 0];[0.6 0.6 0.6],[0 0 0]};
close all
binrange = [0:20:4000];
% binrange = [0:50:580,1200];

bincount = cellfun(@(x) cumsum(histc(x,binrange))'/sum(~isnan(x)),distance_sameprobe,'UniformOutput',false);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 20*0.18]);
hold on;
for i = 1:2
   plot(binrange,cell2mat(bincount(:,i)),'Color',clr{i,1},'LineWidth',0.35);
   plot(binrange,mean(cell2mat(bincount(:,i))),'Color',clr{i,2},'LineWidth',1);
end