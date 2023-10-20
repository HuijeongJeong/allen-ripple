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
avefr = nan(length(cidx),1);
for iS = 1:length(sessionlist)
    iS
    if sum(cidx(session==sessionlist(iS)))==0
        continue;
    end
    load([sdir(sessionlist(iS)),'_cellTable.mat'],'T');
    [in,idx] = ismember(unitid,T.unit_id);
    avefr(in) = T.firing_rate(idx(in));
end

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
% print(fHandle,'-depsc','-painters',[dir,'\depth_distribution.ai']);



%%
layer_n = nan(length(unique(sessionlist)),4,4);
for iS = 1:length(unique(sessionlist))
    for iClst = 1:4
        inclst = cidx==iClst-1;
        if sum(inclst&celltype(:,iCT)&session==sessionlist(iS))==0
           continue;
        end
    layer_n(iS,iClst,:) = histc(layer(inclst&celltype(:,iCT)&session==sessionlist(iS)),binrange_layer);
    end
end
layer_n(isnan(layer_n(:,2,1)),:,:) = [];
%%
close all

fraction_nomod = squeeze(layer_n(:,1,:)./repmat(sum(layer_n(:,1,:),3),1,1,4))*100;
fraction_mod = squeeze(sum(layer_n(:,2:4,:),2)./repmat(sum(sum(layer_n(:,2:4,:),2),3),1,1,4))*100;
inanal = layer_n(:,2:4,:);
inanal = squeeze(sum(inanal,3));
inanal = sum(inanal<5,2)==0;
%%
close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 4]);
hold on;
plot(1:4,fraction_nomod(inanal,:),'Color',[0.6,0.6,0.6],'LineWidth',0.35);
plot(1:4,fraction_mod(inanal,:),'Color',[1,0.6,0.6],'LineWidth',0.35);
errorbar(1:4,mean(fraction_nomod(inanal,:)),std(fraction_nomod(inanal,:))/sqrt(sum(inanal)),...
    'Color','k','CapSize',3);
errorbar(1:4,mean(fraction_mod(inanal,:)),std(fraction_mod(inanal,:))/sqrt(sum(inanal)),...
    'Color','r','CapSize',3);
[h,p,~,stat] = ttest(fraction_nomod(inanal,:),fraction_mod(inanal,:));
p = p*4;
xlim([0.5 4.5])
set(gca,'XTick',1:4,'XTickLabel',{'L2/3','L4','L5','L6'},'XTickLabelRotation',45,...
    'YTick',0:20:80,'XLim',[0 5],'YLim',[0 80],'Box','off','TickDir','out','FontSize',8);
ylabel('% neuron');
dir = 'D:\OneDrive - UCSF\figures\2.allen\revision\figS3';
% print(fHandle,'-depsc','-painters',[dir,'/layer_distribution.ai']);

%%
close all
ct = [cbrewer('qual','Dark2',3);0.6 0.6 0.6];
figure;
hold on;
layer_n(sum(isnan(layer_n(:,2,:)),3)>0,:,:) = nan;
for iL = 1:4
    fraction(iL,:) = squeeze(nansum(layer_n(:,:,iL),1))/sum(squeeze(nansum(layer_n(:,:,iL),1)));
    if iL==1
        l = [2,3];
    else
        l = iL+2;
    end
    avefr_layer(iL) = nanmean(avefr(celltype(:,1) & ismember(layer,l)));
    stdfr_layer(iL) = nanstd(avefr(celltype(:,1) & ismember(layer,l)));
    h = bar(iL,fraction(iL,[2,3,4,1]),'stacked');
    for ic = 1:4
        set(h(ic),'FaceColor',ct(ic,:));
    end
end

