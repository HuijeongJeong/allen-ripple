clearvars; clc; close all;

rng(2)

load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');

layertable = readtable('D:\OneDrive\1.allen-andermann\layer_info.csv');


%% setting up 
k = 3; % number of cluster 

invis = tag.area.vis6;
celltype = [tag.celltype.rs(invis),tag.celltype.fs(invis)];
area = tag.info.structure(invis);
unitid = tag.info.unit_id(invis);
session = tag.info.session_id(invis);
sessionlist = unique(session);
[in,idx] = ismember(unitid,unit_id.vis);
cidx = zeros(sum(invis),1);
cidx(in) = cluster_idx.vis{k-1}(idx(in));

[~,idx] = ismember(unitid,layertable.ecephys_unit_id);
layer = layertable.cortical_layer(idx);
depth = layertable.cortical_depth(idx);

%%
close all
ct = [0.6 0.6 0.6; cbrewer('qual','Dark2',3)];

nIter = 100;
nCell = 100;
    
binrange_depth = 0:0.05:1;
binrange_layer = [2,4,5,6];
iCT = 1;

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 20*0.12]);
hold on;
for iClst = 1:4
    depth_cumsum = NaN(nIter,length(binrange_depth));
    for iIter = 1:nIter
        in = cidx==iClst-1 & celltype(:,iCT);
        depth_sub = depth(randsample(find(in),nCell));
        layer_sub = layer(randsample(find(in),nCell));
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
dir = 'D:\OneDrive - University of California, San Francisco\figures\allen\figS3';
print(fHandle,'-depsc','-painters',[dir,'\depth_distribution.ai']);



%%
layer_n = nan(length(unique(sessionlist)),2,4);
for iS = 1:length(unique(sessionlist))
    for iClst = 1:2
    if iClst==1
       inclst = cidx==0; 
    else
       inclst = cidx>0;
       if sum(inclst&celltype(:,iCT)&session==sessionlist(iS))==0
           continue;
       end
    end
    layer_n(iS,iClst,:) = histc(layer(inclst&celltype(:,iCT)&session==sessionlist(iS)),binrange_layer)/...
        sum(inclst&celltype(:,iCT)&session==sessionlist(iS));
    end
end
layer_n(isnan(layer_n(:,2,1)),:,:) = [];

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 4]);
hold on;
plot(1:4,squeeze(layer_n(:,1,:))*100,'Color',[0.6,0.6,0.6],'LineWidth',0.35);
plot(1:4,squeeze(layer_n(:,2,:))*100,'Color',[1,0.6,0.6],'LineWidth',0.35);
errorbar(1:4,mean(squeeze(layer_n(:,1,:)),1)*100,std(squeeze(layer_n(:,1,:)),[],1)*100/sqrt(size(layer_n,1)),...
    'Color','k','CapSize',3,'Color','k');
errorbar(1:4,mean(squeeze(layer_n(:,2,:)),1)*100,std(squeeze(layer_n(:,2,:)),[],1)*100/sqrt(size(layer_n,1)),...
    'Color','k','CapSize',3,'Color','r');
[h,p,~,stat] = ttest(squeeze(layer_n(:,1,:)),squeeze(layer_n(:,2,:)));
p = p*4;
xlim([0.5 4.5])
set(gca,'XTick',1:4,'XTickLabel',{'L2/3','L4','L5','L6'},'XTickLabelRotation',45,...
    'YTick',0:20:80,'XLim',[0 5],'YLim',[0 80],'Box','off','TickDir','out','FontSize',8);
ylabel('Fraction (%)');
dir = 'D:\OneDrive - University of California, San Francisco\figures\allen\figS3';
print(fHandle,'-depsc','-painters',[dir,'/layer_distribution.ai']);

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
%%

p = [];
pair = {};
for ri = 1:3
    for rj = ri+1:4
        pc = [];
        pairc = {};
        for ci = 1:3
            for cj = ci+1:4
                pc = [pc,chisqNxN(n([ri,rj],[ci,cj]))*36];
                pairc = [pairc,[ri,rj,ci,cj]];
            end
        end
        p = [p;pc];
        pair = [pair;pairc];
    end
end
p(p>1) =1;
h = p<0.05;
%%

p = [];
for i = 1:4
    p = [p;chisqNxN([n(:,i),sum(n(:,[1:i-1,i+1:4]),2)])*4];
end
