clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

win = [-5 5];
binsize = 0.01;
resolution = 10;
mininterval = 0.5;

[cidx,celltype] = deal(cell(nS,1));
iri = cell(nS,6);
% ripplehist = cell(3,4);
data = cell(nS,2,2);
rippletype = {'medial','lateral','global','global_m','global_l'};
numripple = nan(nS,6);
nIter = 10;

for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_ripple_modulated.mat'],'fr_ripple');
    load([sdir(sessionList(iS)),'_cellTable.mat'],'T');
    load([sdir(sessionList(iS)),'_ripples.mat'],'CA1_ripple_classified');
    
    [in,idx] = ismember(T.unit_id,unit_id.vis);
    cidx{iS} = cluster_idx.vis{2}(idx(in));
    [~,idx] = ismember(T.unit_id(in),tag.info.unit_id);
    celltype{iS} = [tag.celltype.rs(idx), tag.celltype.fs(idx)];
    
    totalripple = [CA1_ripple_classified.medial(:,[1,3]);...
        CA1_ripple_classified.lateral(:,[1,3]);
        CA1_ripple_classified.global];
    
    %inter-ripple-interval
    for iRp = 1:3
        iri{iS,iRp} = diff(CA1_ripple_classified.(rippletype{iRp})(:,1));
        if iRp<3
            iri{iS,iRp+3} = diff(sort([CA1_ripple_classified.(rippletype{iRp})(:,1);...
                CA1_ripple_classified.(rippletype{iRp+3})(:,1)]));
            [~,isolatedindex_temp] = isolateRipple(sort([CA1_ripple_classified.(rippletype{iRp})(:,1);...
                CA1_ripple_classified.(rippletype{iRp+3})(:,1)]),mininterval,'on');
            isolatedfraction(iS,iRp) = mean(isolatedindex_temp);
        end
    end
    %%
    iri{iS,6} = diff(sort(totalripple(:,1)));
%     [~,sortidx] = sort(totalripple(:,1));
%     iritemp =  [nan;diff(sort(totalripple(:,1)))];
%     iri{iS,6} = iritemp(
    %%
    % isolate riples
    [~,isolatedindex_total] = isolateRipple(totalripple,mininterval,'on');
    isolatedindex = {isolatedindex_total(1:size(CA1_ripple_classified.medial,1));...
        isolatedindex_total([1:size(CA1_ripple_classified.lateral,1)]+size(CA1_ripple_classified.medial,1))};
    
    % random ripples
    nripple = [size(CA1_ripple_classified.medial,1), size(CA1_ripple_classified.lateral,1)];
%     randomindex = cellfun(@(x) randsapmle(x,
    
    
    % number of total riples & isolated ripples
    numripple(iS,:) = [cellfun(@length,isolatedindex);length(isolatedindex_total);...
        cellfun(@sum,isolatedindex);sum(isolatedindex_total)];
    
%     randomindex = cellfun(@(y) cellfun(@(x) randsample(x,round(x/2)),num2cell(nripple),...
%         'UniformOutput',false),num2cell(1:nIter)','UniformOutput',false);
%     randomindex = cat(1,randomindex{:});

    for i = 1:2
        if ~ismember('hist',fieldnames(fr_ripple))
            spkhist = cellfun(@(y) cell2mat(cellfun(@(x) histc(y,x+[win(1):binsize:win(2)])'*(1/binsize),...
                num2cell(CA1_ripple_classified.(rippletype{i})(:,1)),'UniformOutput',false)),...
                T.spike_time(in),'UniformOutput',false);
            spkhist = cellfun(@(x) x(:,1:end-1),spkhist,'UniformOutput',false);
        else
            spkhist = fr_ripple.hist{i}(in);
        end
        spkave = conv2(cell2mat(cellfun(@(x) nanmean(x),spkhist,'UniformOutput',false)),...
            fspecial('Gaussian',[1 5*resolution],resolution),'same');
        spkave_iso = conv2(cell2mat(cellfun(@(x) nanmean(x(isolatedindex{i},:)),spkhist,'UniformOutput',false)),...
            fspecial('Gaussian',[1 5*resolution],resolution),'same');
        data{iS,i,1} = spkave;
        data{iS,i,2} = spkave_iso;
        
%         for iIter = 1:nIter
%             spkave_random = conv2(cell2mat(cellfun(@(x) nanmean(x(randomindex{iIter,i},:)),spkhist,'UniformOutput',false)),...
%                 fspecial('Gaussian',[1 5*resolution],resolution),'same');
%             data{iS,i,2+iIter} = spkave_random;
%         end
    end
end

%%
intime = fr_ripple.time>=-1.5 & fr_ripple.time<=1.5;
time = fr_ripple.time(intime);

dataz = zscore(cell2mat(cellfun(@(x) x(:,intime),...
        data(:,:,2),'UniformOutput',false)),[],2);
[coef,score,~,~,explained] = pca(dataz);
nPC = find(cumsum(explained)>80,1,'first');

% dataz_original = zscore(cell2mat(cellfun(@(x) x(:,intime),...
%         data(:,:,1),'UniformOutput',false)),[],2);
% [coef_original,score_original] = pca(dataz_original);

% 
% dataz_random = zscore(cell2mat(cellfun(@(x) x(:,intime),...
%         data(:,:,3),'UniformOutput',false)),[],2);
% [coef,score,~,~,explained] = pca(dataz_random);
% nPC = find(cumsum(explained)>80,1,'first');

%%
nC = size(score,1);
probMat = cell(5,1);
for k = 3
kmtmp = NaN(nC,100);
for iIter = 1:100
    kmtmp(:,iIter) = kmeans(score(:,1:nPC),k);
end
probmat = NaN(nC,nC);
for iC = 1:nC
    for jC = iC:nC
        probmat(iC,jC) = sum(kmtmp(iC,:)==kmtmp(jC,:));
        probmat(jC,iC) = sum(kmtmp(iC,:)==kmtmp(jC,:));
    end
end
Z = linkage(probmat,'average','chebychev');
km = cluster(Z,'MaxClust',k);
[~,~,outperm] = dendrogram(Z,0);
probMat{k} = probmat(outperm,outperm);
probIdx = km(outperm);
clusterIdx(:,k) = km;
end
%%
for k= 2:5
   subplot(1,4,k-1);
   imagesc(probMat{k});
   title(['k=',num2str(k)]);
end
colormap('gray')

%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 4]);
k = 3;
c_cluster = cbrewer('qual','Dark2',6);
c_cluster = c_cluster(4:6,:);
axes('Position',axpt(20,20,1,1:19,axpt(1,10,1,2:10,[],[0.02 0.02]),[0.02 0.02]));
imagesc(1,1:nC,probIdx);
ylim([1 nC]);
colormap(gca,c_cluster);
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'YTick',0:1000:3000,'XTick',[]);
ylabel('Neurons','FontSize',8);
axis xy

axes('Position',axpt(20,20,2:20,20,axpt(1,10,1,2:10,[],[0.02 0.02]),[0.02 0.02]));
imagesc(1:nC,1,probIdx');
colormap(gca,c_cluster);
xlim([1 nC]);
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',0:1000:3000,'YTick',[]);
xlabel('Neurons','FontSize',8);


axes('Position',axpt(20,20,2:20,1:19,axpt(1,10,1,2:10,[],[0.02 0.02]),[0.02 0.02]));
imagesc(1:nC,1:nC,probMat{k});
set(gca,'Box','on','XTick',0:1000:3000,'XTickLabel',[],'YTick',0:1000:3000,...
    'YTickLabel',[],'FontSize',8,'LineWidth',0.35);
colormap(gca,'gray');
title(['k=',num2str(k)]);
axis xy
print(fHandle,'-depsc','-painters',...
    'D:\OneDrive - UCSF\figures\2.allen\revision\FigS_stability\custering_result_probmat_isolated.ai');

%%

original_clusteridx = cell2mat(cidx);
for i = 1:3
    for j = 1:3
        prob(i,j) = sum(clusterIdx(:,3)==i&original_clusteridx==j)/sum(original_clusteridx==j);        
    end
end
%%
close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 4]);
hold on;
for i = 1:3
    b = bar(i,prob(:,i),'stacked', 'FaceColor','flat');
    for iclst = 1:3
        b(iclst).CData = c_cluster(iclst,:);
        text(i-0.3,mean([sum(prob(1:iclst-1,i)),sum(prob(1:iclst,i))]),...
            [num2str(round(prob(iclst,i)*1000)/10),'%'],'FontSize',5);
    end
end
set(gca,'XLim',[0.2 3.8],'Box','off','TickDir','out','XTick',1:3,...
    'XTickLabel',{'iAct','dAct','Inh'},'XTickLabelRotation',45,'YTick',0:0.5:1,'YTickLabel',[0,50,100]);
ylabel('% neuron')
print(fHandle,'-depsc','-painters',...
    'D:\OneDrive - UCSF\figures\2.allen\revision\FigS_stability\custering_result_original_isolated.ai');

%%
close all
clr = cbrewer('qual','Dark2',3);

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 12 6]);
intime = fr_ripple.time>=-1.5 & fr_ripple.time<=1.5;
time = fr_ripple.time(intime);
for i = 1:2
    clstidx = cell2mat(cellfun(@(x,y) x(y(:,1)),cidx,celltype,'UniformOutput',false));
    dataz = zscore(cell2mat(cellfun(@(x,y) x(y(:,1),intime),...
        data(:,:,i),repmat(celltype,1,2),'UniformOutput',false)),[],2);
    [~,score] = pca(dataz);
    if i==1
        [~,sortidx] = sortrows([clstidx,...
            max(abs(dataz(:,[100:200,400:500])),[],2)],{'descend','ascend'});
    end
    for iRp = 1:2
        %         axes('Position',axpt(7,1,[1:3]+(iRp-1)*3,1,axpt(2,3,i,1:2)));
        axes('Position',axpt(2,1,iRp,1,axpt(15,1,[1:7]+(i-1)*7,1)));
        hold on;
        imagesc(time,1:length(sortidx),dataz(sortidx,[1:300]+(iRp-1)*300));
        plot([0 0],[0.5 length(sortidx)+0.5],'k:');
        axis xy
        axis tight
        set(gca,'CLim',[-2 2],'XTick',-1:1,'YTick',0:400:1600,'Box','off',...
            'TickDir','out','FontSize',8)
        xlabel('Time (s)');
        if iRp==2
            set(gca,'YTickLabel',[]);
        else
            if i==1
                ylabel('Cortical neuron #');
                title('All ripples');
            else
                set(gca,'YTickLabel',[]);
                title('Isolated ripples');
            end
        end
    end
    
    
end
h = axes('Position',axpt(15,1,15,1));
imagesc(clstidx(sortidx));
axis xy
colormap(h,clr);
set(gca,'Box','off','TickDir','out','XTick',[],'YTick',0:400:1600,'YTickLabel',[]);
print(fHandle,'-depsc','-painters',...
    'D:\OneDrive - UCSF\figures\2.allen\revision\FigS_stability\psth_isolated.ai')

%%
close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 4]);
for iRp = 4:5
    subplot(1,2,iRp-3);
    hold on;
    irihist = cell2mat(cellfun(@(x) histc(x,0:0.01:1.5)',iri(:,iRp),'UniformOutput',false));
    irihist(:,end) = cellfun(@(x) sum(x>=1.5),iri(:,iRp));
    plot(0.005:0.01:1.505,cumsum(irihist,2)./repmat(sum(irihist,2),1,size(irihist,2)),'Color',[0.8 0.8 0.8]);
    plot(0.005:0.01:1.505,mean(cumsum(irihist,2)./repmat(sum(irihist,2),1,size(irihist,2))),'k');
    xlim([0 1.505])
    ylim([0 1]);
    title(titleList{iRp-3});
    if iRp==4
    ylabel('Cumulative fraction (%)');
    xlabel('Inter ripple interval (s)');
    end
    set(gca,'Box','off','TickDir','out','FontSize',8,'XTick',0:0.5:1.5,'YTick',0:0.5:1,'YTickLabel',0:50:100);

end
print(fHandle,'-depsc','-painters',...
    'D:\OneDrive - UCSF\figures\2.allen\revision\FigS_stability\ripple_interval_zoomin.ai');
%%
close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 4]);
iRp = 6;
hold on;
irihist = cell2mat(cellfun(@(x) histc(x,0:0.01:1.5)',iri(:,iRp),'UniformOutput',false));
irihist(:,end) = cellfun(@(x) sum(x>=1.5),iri(:,iRp));
plot(0.005:0.01:1.505,cumsum(irihist,2)./repmat(sum(irihist,2),1,size(irihist,2)),'Color',[0.8 0.8 0.8]);
plot(0.005:0.01:1.505,mean(cumsum(irihist,2)./repmat(sum(irihist,2),1,size(irihist,2))),'k');
xlim([0 1.505])
ylim([0 1]);
title('All ripples');
ylabel('Cumulative fraction (%)');
xlabel('Inter ripple interval (s)');
set(gca,'Box','off','TickDir','out','FontSize',8,'XTick',0:0.5:1.5,'YTick',0:0.5:1,'YTickLabel',0:50:100);
print(fHandle,'-depsc','-painters',...
    'D:\OneDrive - UCSF\figures\2.allen\revision\FigS_stability\ripple_interval_all.ai');

%%
clc
mean(numripple(:,4:6)./numripple(:,1:3))
std(numripple(:,4:6)./numripple(:,1:3))*1.96


% print(fHandle,'-depsc','-painters',...
%     'D:\OneDrive - UCSF\figures\2.allen\revision\FigS_stability\ripple_interval_zoomin.ai');
