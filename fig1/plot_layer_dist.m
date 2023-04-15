clearvars; clc; close all;

rng(2)

load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');

layertable = readtable('D:\OneDrive\1.allen-andermann\layer_info.csv');


%% setting up 
k = 3; % number of cluster 

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

invis = tag.area.vis;
celltype = [tag.celltype.rs(invis),tag.celltype.fs(invis)];
area = tag.info.structure(invis);
unitid = tag.info.unit_id(invis);
[in,idx] = ismember(unitid,unit_id.vis);
cidx = zeros(sum(invis),1);
cidx(in) = cluster_idx.vis{k-1}(idx(in));
session_id = tag.info.session_id(invis);

[~,idx] = ismember(unitid,layertable.ecephys_unit_id);
layer = layertable.cortical_layer(idx);
depth = layertable.cortical_depth(idx);

%%
close all
ct = [0.6 0.6 0.6; cbrewer('qual','Dark2',3)];

nIter = 100;
nCell = 100;
layer_dist = NaN(nIter,4,4);
    
binrange_depth = 0:0.05:1;
binrange_layer = [2,4,5,6];
iCT = 1;

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 20*0.12]);
hold on;
for iClst = 1:4
    in = find(cidx==iClst-1 & celltype(:,iCT) & ismember(session_id,sessionList));
    depth_cumsum = NaN(nIter,length(binrange_depth));
    for iIter = 1:nIter
        depth_sub = depth(randsample(in,nCell));
        layer_sub = layer(randsample(in,nCell));
        depth_cumsum(iIter,:) = cumsum(histc(depth_sub,binrange_depth))/nCell;
        layer_dist(iIter,:,iClst) = histc(layer_sub,binrange_layer)/nCell;
    end
    m_depth = mean(depth_cumsum*100);
    depth_cumsum = sort(depth_cumsum*100);
    s_depth = depth_cumsum([0.05 0.95]*nIter,:);
    fill([binrange_depth flip(binrange_depth)],...
        [s_depth(1,:) flip(s_depth(2,:))],ct(iClst,:),'EdgeColor','none');
    plot(binrange_depth,m_depth,'Color',ct(iClst,:));
end
set(gca,'XTick',0:0.2:1,'YTick',0:20:100,'XLim',[0 1],'YLim',[0 100],'Box','off','TickDir','out','FontSize',8);
ylabel('Cumulative fraction (%)');
xlabel('Normalized depth');
alpha(0.2);
dir = 'D:\OneDrive - UCSF\figures\2.allen\figS3';
print(fHandle,'-depsc','-painters',[dir,'\depth_distribution.ai']);

%%
depthdist = [];
for iClst = 1:4
    depthdist(:,iClst) = histc(depth(cidx==iClst-1 & celltype(:,iCT) & ismember(session_id,sessionList)),...
        binrange_depth);
    for iLayer = 1:4
        n(iClst,iLayer) = sum(cidx==iClst-1 & celltype(:,iCT) & layer==binrange_layer(iLayer) & ismember(session_id,sessionList));
    end
end
depthdist = cumsum(depthdist(1:end-1,:));
% depthdist = cumsum(depthdist(1:end-1,:),1)./repmat(sum(depthdist(1:end-1,:),1),length(binrange_depth)-1,1);
%%
p = [];
k = [];
for i = 1:3
    for j = i+1:4
        [~,ptemp,stat] = kstest2(depth(cidx==i-1 & celltype(:,iCT)& ismember(session_id,sessionList)),...
            depth(cidx==j-1 & celltype(:,iCT)& ismember(session_id,sessionList)));
        p = [p;ptemp*6];
        k = [k;stat];
    end
end
%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 4]);
hold on;
for iClst = 1:4
    data = squeeze(layer_dist(:,:,iClst));
    m_layer = mean(data*100);
    data = sort(data*100);
    s_layer = [data([0.05 0.95]*nIter,:); NaN(1,4)];
    x = [repmat(1:4,2,1);NaN(1,4)];
    plot([1:4]+0.1*(iClst-1),m_layer,'Color',ct(iClst,:));
    plot(x(:)+0.1*(iClst-1),s_layer(:),'Color',ct(iClst,:));
end
set(gca,'XTick',1:4,'XTickLabel',{'L2/3','L4','L5','L6'},'XTickLabelRotation',45,...
    'YTick',0:20:60,'XLim',[0 5],'YLim',[0 60],'Box','off','TickDir','out','FontSize',8);
ylabel('Fraction (%)');
print(fHandle,'-depsc','-painters',[dir,'/layer_distribution.ai']);

