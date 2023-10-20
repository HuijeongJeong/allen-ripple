clearvars; clc; close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric','units');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.global_ripple_number>0 & session_metric.distance_bw_ml_probes>2000);
nS = length(sessionList);

overlap_fraction = nan(nS,4);
data = cell(nS,4);
unitid = cell(nS,1);

win = [-5,5];
analwin = [-1.5,1.5];
bin = 0.01;
resolution = 10;
time = win(1)+bin/2:bin:win(2)-bin/2;
intime = time>=analwin(1) & time<=analwin(2);

for iS = 1:nS
    iS
    clear pupil_data
    load([sdir(sessionList(iS)),'_Events.mat'],'running_speed','pupil_data');
    load([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_win',...
        'spontaneous_anal_win','CA1_ripple_classified');
    load([sdir(sessionList(iS)),'_cellTable.mat'],'T');
    load([sdir(sessionList(iS)),'_ripple_modulated.mat'],'fr_ripple');
    
    %%
    medialidx = find(strcmp(spontaneous_CA1_ripple.relative_location,'medial'));
    lateralidx = find(strcmp(spontaneous_CA1_ripple.relative_location,'lateral'));
    ccf = cell2mat(spontaneous_CA1_ripple.ccf_coordinate);        
    
    %%
    ripples = CA1_ripple_classified.medial
    [pairid,noise_corr_sp,time_sp,pupil_sp(iS),nbin_sp(iS)] =...
        noisecorr_nostim(T.spike_time(invis),spontaneous_win,binsize,unitid_tmp,ripples(:,1),...
        iState,pupil_data.time,pupil_size,running_speed.time,running_speed.velocity_conv,speedlimit);
    
end
%%

datatotal = cellfun(@(x,y) [x,y],data(:,1),data(:,2),'UniformOutput',false);
datatotal = cat(1,datatotal{:});
dataz = zscore(datatotal,[],2);
[~,score,~,~,explained] = pca(dataz);
nPC = find(cumsum(explained)>80,1,'first');
nCell = size(datatotal,1);

%%
nIter = 100;
kmtmp = nan(4,nCell,nIter);
[probMat,probIdx] = deal(cell(4,1));
clusterIdx = nan(nCell,4);
tot_withinss = nan(nIter,5);
for k = 2:5
    k
for iIter = 1:nIter
    [kmtmp(k-1,:,iIter),~,~,D] = kmeans(score(:,1:nPC),k);
    dist = min(D,[],2);
    tot_withinss(iIter,k) = sum(dist);
end

if k <6
probmat = NaN(nCell,nCell);
for iC = 1:nCell
    for jC = iC:nCell
        probmat(iC,jC) = sum(kmtmp(k-1,iC,:)==kmtmp(k-1,jC,:));
        probmat(jC,iC) = sum(kmtmp(k-1,iC,:)==kmtmp(k-1,jC,:));
    end
end

Z = linkage(probmat,'average','chebychev');
km = cluster(Z,'MaxClust',k);
[~,~,outperm] = dendrogram(Z,0);
probMat{k} = probmat(outperm,outperm);
probIdx{k} = km(outperm);
clusterIdx(:,k) = km;
end
end
totss = sum(pdist2(score(:,1:nPC),nanmean(score(:,1:nPC))).^2);
%%
expvar = (totss-tot_withinss)/totss*100;
%%
c_cluster = cbrewer('qual','Dark2',8);
c_cluster(1:3,:) = [];
close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 15 4]);
for k= 2:5

    axes('Position',axpt(20,20,1,1:19,axpt(4,10,k-1,2:10,[],[0.02 0.02]),[0.02 0.02]));
    imagesc(1,1:nCell,probIdx{k});
    ylim([1 nCell]);
    colormap(gca,c_cluster(1:k,:));
    set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
        'YTick',0:1000:3000,'XTick',[]);
    if k==2
        ylabel('Neurons','FontSize',8);
    else
        set(gca,'YTickLabel',[]);
    end
    axis xy
    
    axes('Position',axpt(20,20,2:20,20,axpt(4,10,k-1,2:10,[],[0.02 0.02]),[0.02 0.02]));
    imagesc(1:nCell,1,probIdx{k}');
    colormap(gca,c_cluster(1:k,:));
    xlim([1 nCell]);
    set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
        'XTick',0:1000:3000,'YTick',[]);
    xlabel('Neurons','FontSize',8);
    
    axes('Position',axpt(20,20,2:20,1:19,axpt(4,10,k-1,2:10,[],[0.02 0.02]),[0.02 0.02]));
    imagesc(1:nCell,1:nCell,probMat{k});
    set(gca,'Box','on','XTick',0:1000:3000,'XTickLabel',[],'YTick',0:1000:3000,...
        'YTickLabel',[],'FontSize',8,'LineWidth',0.35);
    colormap(gca,'gray');
    title(['k=',num2str(k)]);
    axis xy
end
print(fHandle,'-dtiff','-r600','D:\OneDrive - UCSF\data\allen\revision\clustering_transverse.tif');
%%
[~,idx] = ismember(cell2mat(unitid),unit_id.vis);
clusterIdx_original = cluster_idx.vis{2}(idx);

close all
for k = 2:3
    subplot(1,2,k-1);
probcluster = nan(k,3);
for i = 1:k
    for j = 1:3
        probcluster(i,j) = sum(clusterIdx(:,k)==i & clusterIdx_original==j)/sum(clusterIdx_original==j);
    end
end
imagesc(probcluster);
end
%%
%%
close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 7 4]);
for k = 2:3
    axes('Position',axpt(2,1,k-1,1,axpt(10,10,2:10,1:9),[0.05 0.05]));
    hold on;
    probcluster = nan(k,3);
    for i = 1:k
        for j = 1:3
            probcluster(i,j) = sum(clusterIdx(:,k)==i & clusterIdx_original==j)/sum(clusterIdx_original==j);
        end
    end
    
    for i = 1:3
        b = bar(i,probcluster(:,i),'stacked', 'FaceColor','flat');
        for iclst = 1:k
            b(iclst).CData = c_cluster(iclst,:);
            text(i-0.3,mean([sum(probcluster(1:iclst-1,i)),sum(probcluster(1:iclst,i))]),...
                [num2str(round(probcluster(iclst,i)*1000)/10),'%'],'FontSize',5);
        end
    end
    set(gca,'XLim',[0.2 3.8],'Box','off','TickDir','out','XTick',1:3,...
        'XTickLabel',{'iAct','dAct','Inh'},'XTickLabelRotation',45,'YTick',0:0.5:1,'YTickLabel',[0,50,100]);
    if k==2
    ylabel('% neuron')
    else
       set(gca,'YTickLabel',[]); 
    end
end
print(fHandle,'-dtiff','-r600','D:\OneDrive - UCSF\data\allen\revision\clustering_result_original_transverse.tif');