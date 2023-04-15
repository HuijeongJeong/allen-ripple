clearvars; clc; close all;

rng(2)

load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_vis.mat');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.global_ripple_number>0 & session_metric.distance_bw_ml_probes>2000);
nS = length(sessionList);

nIter = 10;
nK = 20;

refList = {'medial';'lateral'};
    
ylimit = [-2 2; -2 3];

intime = timespk>=-1.5 & timespk<=1.5;
[data,celltype,modref,unitid] = deal(cell(2,1));
for iType = 1:2
    modnew = cellfun(@(x,y,z) x|y,mod(:,iType,1),mod(:,iType,2),'UniformOutput',false);
    data{iType} = cellfun(@(x,y,z) [x(z,intime),y(z,intime)],spkconv(:,iType,1),...
        spkconv(:,iType,2),modnew,'UniformOutput',false);
    celltype{iType} = cellfun(@(x) ones(size(x,1),1)*iType,data{iType},'UniformOutput',false);
    modref{iType} = cellfun(@(x,y,z) [x(z),y(z)],mod(:,iType,1),mod(:,iType,2),modnew,'UniformOutput',false);
    unitid{iType} = cellfun(@(x,y) x(y),unitID(:,iType),modnew,'UniformOutput',false);
end
unitid = cell2mat(cat(1,unitid{:}));
modref = cell2mat(cat(1,modref{:}));
datatotal = cell2mat(cat(1,data{:}));
celltype = cell2mat(cat(1,celltype{:}));
dataz = zscore(datatotal,[],2);
[coef,score,~,~,explained] = pca(dataz);
nPC = find(cumsum(explained)>80,1,'first');

%%
dir = 'D:\OneDrive - University of California, San Francisco\figures\allen\figS2';

clr = {'r';'b'};
time = timespk(intime);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 12 3.5]);
for iPC = 1:3
    axes('Position',axpt(3,1,iPC,1,axpt(1,10,1,1:9)))
    plot(time,coef(1:length(time),iPC),'r');
    hold on;
    plot(time,coef(length(time)+1:end,iPC),'b');
    plot([0 0],[-1.5 1.5],'k:');
    plot([-1.5 1.5],[0 0],'Color',[0.8 0.8 0.8]);
    ylim([-0.1 0.13])
    xlim([-1.5 1.5]);
    set(gca,'XTick',-1.5:1.5:1.5,'YTick',-0.1:0.1:0.1,'Box','off','TickDir','out',...
        'FontSize',8,'LineWidth',0.35);
    if iPC==1
        xlabel('Time from ripple onset (s)');
        ylabel('Activity (a.u.)');
    else
        set(gca,'YTickLabel',[]);
    end
    title(['PC ',num2str(iPC), ' (',num2str(round(explained(iPC)*10)/10),'%)']);
end
print(fHandle,'-depsc','-painters',[dir,'\cluster_PC.ai']);


%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3.5 3.5]);
plot(0:length(explained),cumsum([0;explained]),'k','LineWidth',0.5);
hold on;
plot([9 9],[0 100],'k--','LineWidth',0.35);
plot([0 100],[80 80],'k:','LineWidth',0.35)
scatter(9,sum(explained(1:9)),5,'r','filled')
xlim([0 45])
xlabel('Number of PC');
ylabel('Explained variance (%)');
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',[0 9 15 30 45],'YTick',0:20:100);
print(fHandle,'-depsc','-painters',[dir,'\cluster_PC_explainedvariance.ai']);

%%
[tot_withinss,clst_density,sil,sil_pca,gap_stat] = deal(NaN(nIter,nK));
for iIter = 1:100
    iIter
    a = NaN(nK,1);
    for k = 1:nK
        [km,~,~,D] = kmeans(score(:,1:nPC),k);
        
        % for silhouette value
        sil(iIter,k) = nanmean(silhouette(dataz,km));
        sil_pca(iIter,k) = nanmean(silhouette(score(:,1:nPC),km));
        dist = min(D,[],2);
        
        % for elbow method
        tot_withinss(iIter,k) = sum(dist); % within cluster sum of squares
        
        % for clustering density
%         if k==1
%             clst_density(iIter,k) = 1;
%         else
%             if k==2
%                 a(k) = 1-3/(4*size(dataz,2));
%             else
%                 a(k) = a(k-1)+(1-a(k-1))/6;
%             end
%             clst_density(iIter,k) = tot_withinss(iIter,k)/(a(k)*tot_withinss(iIter,k-1));
%         end
    end
end
sil(:,1) = 0;
totss = sum(pdist2(score(:,1:nPC),nanmean(score(:,1:nPC))).^2);
expvar = (totss-tot_withinss)/totss*100;

%%
close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 4]);
errorbar(1:20,nanmean(expvar),nanstd(expvar)/sqrt(nIter),'k','CapSize',2);
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',0:5:20,...
    'YTick',0:20:90,'XLim',[0 21],'YLim',[-5 85]);
xlabel('Number of clusters (k)');
ylabel('Explained variance (%)');
print(fHandle,'-depsc','-painters',[dir,'\explainedvariance_kcluster.ai']);

%%
[probMat,clusterIdx,probIdx] = deal(cell(4,1));
nC = size(score,1);
for k = 2:5
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
    probMat{k-1} = probmat(outperm,outperm);
    probIdx{k-1} = km(outperm);
    clusterIdx{k-1} = km;
end
%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 15 4]);
for k = 2:5
    if k<3
        c_cluster = cbrewer('qual','Dark2',3);
    else
    c_cluster = cbrewer('qual','Dark2',k);
    end
    axes('Position',axpt(20,20,1,1:19,axpt(4,10,k-1,2:10,[],[0.02 0.02]),[0.02 0.02]));
    imagesc(1,1:nC,probIdx{k-1});
    ylim([1 nC]);
    colormap(gca,c_cluster);
    set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
        'YTick',0:1000:3000,'XTick',[]);
    if k==2
        ylabel('Neurons','FontSize',8);
    else
        set(gca,'YTickLabel',[]);
    end
    axis xy
    
    axes('Position',axpt(20,20,2:20,20,axpt(4,10,k-1,2:10,[],[0.02 0.02]),[0.02 0.02]));
    imagesc(1:nC,1,probIdx{k-1}');
    colormap(gca,c_cluster);
    xlim([1 nC]);
    set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
        'XTick',0:1000:3000,'YTick',[]);
    xlabel('Neurons','FontSize',8);
    
    
    axes('Position',axpt(20,20,2:20,1:19,axpt(4,10,k-1,2:10,[],[0.02 0.02]),[0.02 0.02]));
    imagesc(1:nC,1:nC,probMat{k-1});
    set(gca,'Box','on','XTick',0:1000:3000,'XTickLabel',[],'YTick',0:1000:3000,...
        'YTickLabel',[],'FontSize',8,'LineWidth',0.35);
    colormap(gca,'gray');
    title(['k=',num2str(k)]);
    axis xy
end
print(fHandle,'-depsc','-painters',[dir,'\prob_matrix_kcluster.ai']);
%%
close all;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 3]);
set(gca,'CLim',[0 100]);
colormap(gca,'gray');
h = colorbar;
h.XTick = [0,50,100];
h.Box = 'off';
axis off
print(fHandle,'-depsc','-painters',[dir,'\prob_matrix_kcluster_colorbar.ai']);


% 
% unit_id.vis = unitid;
% cluster_idx.vis = clusterIdx;
% cluster_k.vis = 2:5;
% proMat.vis = probMat;
% 
% save('D:\heejeong\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat',...
%     'unit_id','cluster_idx','cluster_k','probMat');
% 

