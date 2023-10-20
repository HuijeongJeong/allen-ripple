clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

iState = 1;
binsize = 2;
speedlimit = 2;

k = 3;
[noisecorr_sp,withinarea] = deal(cell(nS,k*(k+1)/2+k+1));
[time_sp,pupil_sp,nbin_sp] = deal(NaN(nS,2));
clusteridx = cell(nS,1);

for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_cellTable.mat'],'T');
    load([sdir(sessionList(iS)),'_ripples.mat'],'CA1_ripple_classified',...
        'spontaneous_win','spontaneous_anal_win');
    load([sdir(sessionList(iS)),'_Events.mat'],'running_speed','pupil_data');
    
    invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis & tag.celltype.rs));
    nCell = sum(invis);
    unitid_tmp = T.unit_id(invis);
    [~,idx] = ismember(unitid_tmp,tag.info.unit_id);
    ccf = tag.info.ccf(idx,:);
    structure = tag.info.structure(idx);
    
    pupil_size = pupil_data.pupil_height/2.*pupil_data.pupil_width/2*pi;
    noise_corr_sp = NaN(nCell*(nCell-1)/2,2);
    
    rippletime = sort([CA1_ripple_classified.medial(:,1);...
        CA1_ripple_classified.lateral(:,1); CA1_ripple_classified.global(:,1)]);
    %%
    window = [spontaneous_win(1),rippletime(round(length(rippletime)/4));...
        rippletime(end-round(length(rippletime)/4)+1),spontaneous_win(2)];
    %%
    for iwin = 1:2
        [pairid,noise_corr_sp(:,iwin),time_sp(iS,iwin),pupil_sp(iS,iwin),nbin_sp(iS,iwin)] =...
            noisecorr_nostim(T.spike_time(invis),window(iwin,:),binsize,unitid_tmp,rippletime,...
            iState,pupil_data.time,pupil_size,running_speed.time,running_speed.velocity_conv,speedlimit);
    end
    
    %%
    ccfout = sum(ccf<0,2)>0;
    ccf(ccfout,:) = NaN;
    pair_distance = pdist(ccf)';
    
    
    [in,idx] = ismember(unitid_tmp,unit_id.vis);
    clusteridx{iS} = zeros(length(idx),1);
    clusteridx{iS}(in) = cluster_idx.vis{k-1}(idx(in));
    
    [~,idx] = ismember(pairid,unitid_tmp);
    cidxtmp = clusteridx{iS}(idx);
    areatmp = structure(idx);
    kClst = 1;
    for iClst = 1:k+1
        for jClst = iClst:k+1
            in = cidxtmp(:,1)==iClst-1 & cidxtmp(:,2)==jClst-1;
            noisecorr_sp{iS,kClst} = noise_corr_sp(in,:);
            withinarea{iS,kClst} = cellfun(@(x,y) strcmp(x,y),areatmp(in,1),areatmp(in,2));
            kClst = kClst+1;
        end
    end
end

%%
nmod = cell2mat(cellfun(@(x) [sum(x==1),sum(x==2),sum(x==3)],clusteridx,'UniformOutput',false));
out = sum(nmod==0,2)>0;
x = [1 1 1 1 2 2 2 3 3 4];
y = [1 2 3 4 2 3 4 3 4 4];
clstlist = {'Nomod';'iAct';'dAct';'Inh'};
data = nan(nS,length(x),2);
c = cbrewer('div','RdBu',50);
c([1:5, 46:50],:) = [];

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.35 20/6]);
for iwin = 1:2
    data(~out,:,iwin) = cellfun(@(x) nanmean(x(:,iwin)),noisecorr_sp(~out,:));
    msquare = nan(4,4,nS);
    m = squeeze(data(:,:,iwin));
    for i = 1:10
       msquare(x(i),y(i),:) = m(:,i); 
       msquare(y(i),x(i),:) = m(:,i); 
    end
    h(iwin) = axes('Position',axpt(2,1,iwin,1,axpt(1,10,1,1:9)));
    imagesc(squeeze(nanmean(msquare,3)));
    axis xy;
    colormap(c)
    set(gca,'CLim',[-0.15 0.15],'XTick',1:4,'XTickLabel',clstlist,'YTick',1:4,...
        'YTickLabel',clstlist,'Box','off','TickDir','out','XTickLabelRotation',45,...
        'FontSize',7);
    if iwin==1
        title('First quadrant')
    else
        title('Last quadrant');
        set(gca,'YTickLabel',[]);
    end
end

[~,p] = ttest(data(:,:,1),data(:,:,2));
for i = 1:10
    if p(i)<0.05
       text(h(1),x(i),y(i),'*');
       text(h(1),y(i),x(i),'*');
       text(h(2),x(i),y(i),'*');
       text(h(2),y(i),x(i),'*');
    end
end
print(fHandle,'-depsc','-painters',...
    'D:\OneDrive - UCSF\figures\2.allen\revision\FigS_stability\spontaneous_corr_firstlast.ai');

