clearvars; clc; close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
layertable = readtable('D:\OneDrive\1.allen-andermann\layer_info.csv');
sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

ripplelist = {'medial','lateral','global'};

winbase = [-2 -1];
winripple = [-0.5 0.5];

[clstidx,layer,depth,selectivity,pselec] = deal(cell(nS,1));
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_Events.mat'],'natural_movie');
    load([sdir(sessionList(iS)),'_ripples.mat'],'CA1_ripple_classified');
    
    invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis & tag.celltype.rs));
    [clstidx{iS},layer{iS},depth{iS}] = deal(zeros(sum(invis),1));

    [in,idx] = ismember(T.unit_id(invis),unit_id.vis);
    clstidx{iS}(in) = cluster_idx.vis{2}(idx(in));
    
    [in,idx] = ismember(T.unit_id(invis),layertable.ecephys_unit_id);
    layer{iS}(in) = layertable.cortical_layer(idx(in));
    depth{iS}(in) = layertable.cortical_depth(idx(in));

    [selectivity{iS},pselec{iS}] = deal(nan(sum(invis),3));
    for iRp = 1:3
        spkbase = cellfun(@(y) cellfun(@(x) sum(y>=x+winbase(1) & y<=x+winbase(2)),...
            num2cell(CA1_ripple_classified.(ripplelist{iRp})(:,1))),...
            T.spike_time(invis),'UniformOutput',false);
        spkripple = cellfun(@(y) cellfun(@(x) sum(y>=x+winripple(1) & y<=x+winripple(2)),...
            num2cell(CA1_ripple_classified.(ripplelist{iRp})(:,1))),...
            T.spike_time(invis),'UniformOutput',false);
        selectivity{iS}(:,iRp) = cellfun(@(x,y)(mean(x)-mean(y))/sqrt(std(x)^2+std(y)^2),spkripple,spkbase);
        [~,pselec{iS}(:,iRp)] = cellfun(@(x,y) ttest(x,y),spkripple,spkbase);
    end    
end

%%
% close all
c_cluster = [1 0.6 0.6; 0.6 0.6 1; 0.6 0.6 0.6];

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4.5 4]);
iRp = 3;
hold on;
for iL = 1:4
    if iL==1
        inlayer = cellfun(@(x) ismember(x,[2,3]),layer,'UniformOutput',false);
    else
        inlayer = cellfun(@(x) x==iL+2,layer,'UniformOutput',false);
    end
    
    num = sum(cell2mat(cellfun(@(x,y,z) sum([x(y,iRp)<0.05 & z(y,iRp)>0,...
        x(y,iRp)<0.05 & z(y,iRp)<0,x(y,iRp)>0.05]),pselec,inlayer,selectivity,'UniformOutput',false)));
    fraction = num/sum(num);
   
    b = bar(iL,fraction,'stacked','FaceColor','flat');
    for i = 1:3
        b(i).CData = c_cluster(i,:);
        text(iL-0.3,mean([sum(fraction(1:i-1)),sum(fraction(1:i))]),...
            [num2str(num(i)),newline,'(',num2str(round(fraction(i)*1000)/10),'%)'],'FontSize',5);
    end
end
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',1:4,...
    'XTickLabel',{'L2/3','L4','L5','L6'},'XTickLabelRotation',45,'YTick',0:0.5:1,...
    'YTickLabel',0:50:100,'XLim',[0.25 4.75]);
ylabel('% neurons');
alpha(0.5)
print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\data\allen\revision\ripple_mod_selectivity.ai');

%%
% close all
clr = [cbrewer('qual','Dark2',3);0.6 0.6 0.6];

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4.5 4]);
hold on;
for iL = 1:4
    if iL==1
        inlayer = cellfun(@(x) ismember(x,[2,3]),layer,'UniformOutput',false);
    else
        inlayer = cellfun(@(x) x==iL+2,layer,'UniformOutput',false);
    end
    
    num = sum(cell2mat(cellfun(@(x,y) sum([x(y)==1,x(y)==2,x(y)==3,x(y)==0]),...
        clstidx,inlayer,'UniformOutput',false)));
    fraction = num/sum(num);
   
    b = bar(iL,fraction,'stacked','FaceColor','flat');
    for i = 1:4
        b(i).CData = clr(i,:);
        text(iL-0.3,mean([sum(fraction(1:i-1)),sum(fraction(1:i))]),...
            [num2str(num(i)),newline,'(',num2str(round(fraction(i)*1000)/10),'%)'],'FontSize',5);
    end
end
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',1:4,...
    'XTickLabel',{'L2/3','L4','L5','L6'},'XTickLabelRotation',45,'YTick',0:0.5:1,...
    'YTickLabel',0:50:100,'XLim',[0.25 4.75]);
ylabel('% neurons');
alpha(0.5)
print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\data\allen\revision\ripple_mod_glm.ai');

%%
clr = cbrewer('qual','Dark2',3);
clrdim = cbrewer('qual','Pastel2',3);
fractionlayer = cell(3,1);
out = false(nS,1);
for iClst = 1:3
    fractionlayer{iClst} = nan(nS,4);
    axes('Position',axpt(3,1,iClst,1));
    hold on;
    inclst = cellfun(@(x) x==iClst,clstidx,'UniformOutput',false);
    out = out | cellfun(@sum,inclst)<5;
    for iL = 1:4
        if iL==1
            inlayer = cellfun(@(x) ismember(x,[2,3]),layer,'UniformOutput',false);
        else
            inlayer = cellfun(@(x) x==iL+2,layer,'UniformOutput',false);
        end
        fractionlayer{iClst}(:,iL) = cellfun(@(x,y) sum(x&y)/sum(x),inclst,inlayer);        
    end
end

close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 4]);
hold on;
for iClst = 1:3
%     plot(1:4,fractionlayer{iClst}(~out,:),'Color',clrdim(iClst,:));
    errorbar([1:4]+(iClst-1)*0.1,mean(fractionlayer{iClst}(~out,:)),...
        std(fractionlayer{iClst}(~out,:))/sqrt(sum(~out)),'Color',clr(iClst,:),'CapSize',3);
end
xlim([0.1 5.1])
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',1:4,...
    'XTickLabel',{'L2/3','L4','L5','L6'},'XTickLabelRotation',45,'YTick',0:0.2:0.6,'YLim',[0,0.6],...
    'YTickLabel',0:20:60);
ylabel('% neurons');
print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\data\allen\revision\clst_layer.ai');

data = cat(3,fractionlayer{:});
data = data(~out,:,:);
data = permute(data,[2,1,3]);
data = data(:,:)';
groupidx = repmat([1,2,3],sum(~out),1);
simple_mixed_anova(data,groupidx(:));