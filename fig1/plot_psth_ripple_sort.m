clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);


win = [-2 2];
binsize = 0.01;
resolution = 10;

[cidx,celltype] = deal(cell(nS,1));
data = cell(nS,3);
data_g = cell(nS,2);

for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_ripple_modulated.mat'],'fr_ripple');
    load([sdir(sessionList(iS)),'_cellTable.mat'],'T');
    load([sdir(sessionList(iS)),'_ripples.mat'],'CA1_ripple_classified');
    
    [in,idx] = ismember(T.unit_id,unit_id.vis);
    data(iS,:) = cellfun(@(x) x(in,:),fr_ripple.conv,'UniformOutput',false);
    cidx{iS} = cluster_idx.vis{2}(idx(in));
    [in,idx] = ismember(T.unit_id(in),tag.info.unit_id);
    celltype{iS} = [tag.celltype.rs(idx(in)), tag.celltype.fs(idx(in))];
        
    spkhist = cellfun(@(y) cell2mat(cellfun(@(x) histc(y,x+[win(1):binsize:win(2)])',...
        num2cell(CA1_ripple_classified.global(:,1)),'UniformOutput',false)),...
        T.spike_time(in),'UniformOutput',false);
    spkhist = cellfun(@(x) x(:,1:end-1),spkhist,'UniformOutput',false);
    dlead = ismember(CA1_ripple_classified.global(:,1),CA1_ripple_classified.global_m(:,1));
    for i = 1:2
        if i==1
            in = dlead;
        else
            in = ~dlead;
        end
        spkave = conv2(cell2mat(cellfun(@(x) nanmean(x(in,:)/binsize),spkhist,'UniformOutput',false)),...
            fspecial('Gaussian',[1 5*resolution],resolution),'same');
        data_g{iS,i} = spkave;
    end
    
    
end

%%
time = fr_ripple.time;
dataplot = cellfun(@(x) x(:,time>=-1.5 & time<=1.5),data(:,1:2),'UniformOutput',false);
m = mean(cell2mat(dataplot),2);
s = std(cell2mat(dataplot),[],2);
dataplotz = zscore(cell2mat(dataplot),[],2);

[coef,score,~,~,explained] = pca(dataplotz);
nPC = find(cumsum(explained)>80,1,'first');
cidx = cell2mat(cidx);
celltype = cell2mat(celltype);

%%
titlelist = {'dCA1 ripple','iCA1 ripple'};
ct = cbrewer('qual','Dark2',3);
timeplot = time(time>=-1.5 & time<=1.5);
d = dataplotz(celltype(:,1),:);
c = cidx(celltype(:,1));
[~,sortidx] = sortrows([cidx(celltype(:,1)),max(abs(d(:,[100:200,400:500])),[],2)],...
    {'descend','ascend'});

%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20/5*2 20/3]);
for i = 1:2
    axes('Position',axpt(13,1,[1:5]+(i-1)*5,1));
    hold on;
    imagesc(timeplot,1:size(d,1),d(sortidx,[1:length(timeplot)]+(i-1)*length(timeplot)));
    plot([0 0],[0.5 size(d,1)+0.5],'k:');
    set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,...
        'CLim',[-2 2],'YLim',[0.5 size(d,1)+0.5],'XLim',[-1.5 1.5],'XTick',-1:1:1,'YTick',0:400:1600);
    if i==2
        set(gca,'YTickLabel',[]);
    else
        ylabel('Neurons');
    end
    title(titlelist{i});
    xlabel('Time (s)');
end

h(1) = axes('Position',axpt(13,1,11:12,1));
s = score(celltype(:,1),1:3);
imagesc(1:3,1:size(d,1),s(sortidx,:));
set(h(1),'YTickLabel',[],'Fontsize',7,'Box','off','TickDir','out');
colormap(h(1),'gray')
xlabel('PC');
axis xy

h(2) = axes('Position',axpt(13,1,13,1));
s = score(celltype(:,1),1:3);
imagesc(1,1:size(d,1),c(sortidx,:));
colormap(h(2),ct);
set(h(2),'YTickLabel',[],'Fontsize',7,'Box','off','TickDir','out','XTick',[]);
axis xy

cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig1')
% print(fHandle,'-depsc','-painters','psth_vis.ai');
%%
dataplotz_gtotal = (cell2mat(data(:,3))-repmat(m,1,length(time)))./repmat(s,1,length(time));
d_gtotal = dataplotz_gtotal(celltype(:,1),:);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20/5 20/3]);
hold on;
imagesc(time,1:size(d,1),d_gtotal(sortidx,:));
plot([0 0],[0.5 size(d,1)+0.5],'k:');
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,...
    'CLim',[-2 2],'YLim',[0.5 size(d,1)+0.5],'XLim',[-1.5 1.5],'XTick',-1:1:1);
ylabel('Neurons');
title('Gloabal');
xlabel('Time (s)');

cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig1')
print(fHandle,'-depsc','-painters','psth_vis_global.ai');

%%
dataplotz_g = cell2mat(data_g);
dataplotz_g = (dataplotz_g-repmat(m,1,size(dataplotz_g,2)))./repmat(s,1,size(dataplotz_g,2));
d_g = dataplotz_g(celltype(:,1),:);
%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20/5*2 20/3]);
for i = 1:2
    axes('Position',axpt(13,1,[1:5]+(i-1)*5,1));
    hold on;
    imagesc(win(1)+binsize/2:binsize:win(2),1:size(d_g,1),d_g(sortidx,[1:diff(win)/binsize]+(i-1)*diff(win)/binsize));
    plot([0 0],[0.5 size(d,1)+0.5],'k:');
    set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,...
        'CLim',[-2 2],'YLim',[0.5 size(d,1)+0.5],'XLim',[-1.5 1.5],'XTick',-1:1:1);
    if i==2
        set(gca,'YTickLabel',[]);
    else
        ylabel('Neurons');
    end
    title(titlelist{i});
    xlabel('Time (s)');
end

h(1) = axes('Position',axpt(13,1,11:12,1));
s = score(celltype(:,1),1:3);
imagesc(1:3,1:size(d,1),s(sortidx,:));
set(h(1),'YTickLabel',[],'Fontsize',7,'Box','off','TickDir','out');
colormap(h(1),'gray')
xlabel('PC');
axis xy

h(2) = axes('Position',axpt(13,1,13,1));
s = score(celltype(:,1),1:3);
imagesc(1,1:size(d,1),c(sortidx,:));
colormap(h(2),ct);
set(h(2),'YTickLabel',[],'Fontsize',7,'Box','off','TickDir','out','XTick',[]);
axis xy

cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig1')
print(fHandle,'-depsc','-painters','psth_vis_global_detail.ai');

%% colorbar
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.15 20*0.18]);
hc = colorbar('North');
% colormap(c)
set(gca,'CLim',[-2 2]);
set(hc,'Box','off','TickDir','out','FontSize',6,'XTick',-2:2);
axis off
cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig1');
print(fHandle,'-depsc','-painters','colorbar.ai')


%%

clr = {'r';'b'};
titleList = {'I act';'D act';'Inh'};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 2.5]);
for iclst = 1:3
    subplot(1,3,iclst);
    hold on;
    for i = 1:2
        dc = d(c==iclst,[1:length(timeplot)]+(i-1)*length(timeplot));
        m = nanmean(dc);
        s = nanstd(dc)/sqrt(sum(c==iclst));
        fill([timeplot,flip(timeplot)],[m+s flip(m-s)],clr{i},'EdgeColor','none');
        plot(timeplot,m,'Color',clr{i});
    end
    dc = d_gtotal(c==iclst,:);
    m = nanmean(dc);
        s = nanstd(dc)/sqrt(sum(c==iclst));
    fill([time,flip(time)],[m+s flip(m-s)],'k','EdgeColor','none');
    plot(time,m,'k');
    
    plot([0 0],[-1.5 2.5],'k:');
    set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,...
        'XTick',-1:1:1,'YTick',-1:1:2,'XLim',[-1.5 1.5],'YLim',[-1.5 2.5]);
    if iclst==1
        ylabel('Firing rate (z-score)');
    else
        set(gca,'YTickLabel',[]);
    end
    xlabel('Time (s)');
    alpha(0.2);
    title(titleList{iclst});
end
print(fHandle,'-dtiff','-r600','ave_psth_vis.tif');

 