clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');

k = 3;
sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);
celltype = [tag.celltype.rs, tag.celltype.fs];

[sigcorr_dg,noisecorr_dg,noisecorr_sp] = deal(cell(nS,k*(k+1)/2+k+1));
[ntrial_dg,time_sp,pupil_dg,pupil_sp,nbin_sp] = deal(NaN(nS,4));
[clusteridx,distance,withinarea,areas] = deal(cell(nS,1));

speedlimit = 2;
binsize = 2;

threshold = [2 5];   %threshold for ripple beginning/end and peak
duration = [30 20 250];    %min inter-ripple inter, min ripple length, max ripple length
Fs = 1250;

%%
typeList = {'rs';'fs'};
iType = 1;
for iS = 1:nS
    iS
    clear pupil_data
    load([sdir(sessionList(iS)),'_Events.mat'],'running_speed','pupil_data');
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_ripples.mat'],'filtered_lfp','spontaneous_anal_win','spontaneous_win');

    % detect ripples
    win = spontaneous_win;
    start = running_speed.time(diff([0;running_speed.immobile])==1);
    stop = running_speed.time(diff(running_speed.immobile)==-1);
    if length(stop)<length(start)
        immobile = [start,[stop;win(2)]];
    else
        immobile = [start,stop];
    end
    ripples = cell(2,1);
    for iProbe = 1:2
        if iscell(filtered_lfp.time{iProbe})
            ripple_tmp = cell2mat(cellfun(@(x,y) FindRipples_HJ([x',y],...
                'thresholds',threshold,'durations',duration,'frequency',Fs),...
                filtered_lfp.time{iProbe},filtered_lfp.lfp{iProbe},'UniformOutput',false));
        else
            ripple_tmp = FindRipples_HJ([filtered_lfp.time{iProbe}',filtered_lfp.lfp{iProbe}],...
                'thresholds',threshold,'durations',duration,'frequency',Fs);
        end
        inImmobile = logical(sum(cell2mat(cellfun(@(x) ripple_tmp(:,1)>x(1) & ripple_tmp(:,3)<x(2),...
            mat2cell(immobile,ones(size(immobile,1),1),2),'UniformOutput',false)'),2));
        ripples{iProbe} = ripple_tmp(inImmobile,:);
    end

    %%

    invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis & celltype(:,iType)));
    unitid_tmp = T.unit_id(invis);
    [~,idx] = ismember(unitid_tmp,tag.info.unit_id);
    structure = tag.info.structure(idx);
    ccf = tag.info.ccf(idx,:);

    [in,idx] = ismember(unitid_tmp,unit_id.vis);
    clusteridx{iS} = zeros(length(idx),1);
    clusteridx{iS}(in) = cluster_idx.vis{k-1}(idx(in));

    %%

    nCell = sum(invis);
    pupil_size = pupil_data.pupil_height/2.*pupil_data.pupil_width/2*pi;
    noise_corr_sp = NaN(nCell*(nCell-1)/2,4);
    ripples = cell2mat(ripples);
    for iState = 1:4
        [pairid,noise_corr_sp(:,iState),time_sp(iS,iState),pupil_sp(iS,iState),nbin_sp(iS,iState)] =...
            noisecorr_nostim(T.spike_time(invis),spontaneous_win,binsize,unitid_tmp,ripples(:,1),...
            iState,pupil_data.time,pupil_size,running_speed.time,running_speed.velocity_conv,speedlimit);
    end
    ccfout = sum(ccf<0,2)>0;
    ccf(ccfout,:) = NaN;
    pair_distance = pdist(ccf)';

    [~,idx] = ismember(pairid,unitid_tmp);
    cidxtmp = clusteridx{iS}(idx);
    areatmp = structure(idx);
    kClst = 1;
    for iClst = 1:k+1
        for jClst = iClst:k+1
            in = cidxtmp(:,1)==iClst-1 & cidxtmp(:,2)==jClst-1;
            noisecorr_sp{iS,kClst} = noise_corr_sp(in,:);
            distance{iS,kClst} = pair_distance(in);
            withinarea{iS,kClst} = cellfun(@(x,y) strcmp(x,y),areatmp(in,1),areatmp(in,2));
            areas{iS,kClst} = categorical(areatmp(in,:)); 
            kClst = kClst+1;
        end
    end
end

%%
nmod = cell2mat(cellfun(@(x) [sum(x==1),sum(x==2),sum(x==3)],clusteridx,'UniformOutput',false));
out = sum(nmod==0,2)>0;

x = [1 1 1 1 2 2 2 3 3 4];
y = [1 2 3 4 2 3 4 3 4 4];

stateList = {'all';'Immobile-low pupil';'Immobile-high pupil';'Mobile'};
dataList = {'dg_signal';'dg_noise';'sp_noise'};
c = cbrewer('div','RdBu',50);
c([1:5, 46:50],:) = [];
%% colorbar
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.15 20*0.18]);
hc = colorbar('North');
colormap(c)
set(gca,'CLim',[-0.15 0.15]);
set(hc,'Box','off','TickDir','out','FontSize',6,'XTick',-0.15:0.15:0.15);
axis off
% cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig2');
% print(fHandle,'-depsc','-painters','colorbar.ai')

%%
iData = 3;
data = NaN(nS,length(x),4);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20/2 20/6]);
for iState = 2:4
    switch iData
        case 1
            in = sum(ntrial_dg(:,2:4)>30,2)==3 & ~out;
            data(in,:,iState) = cellfun(@(x) nanmean(x(:,iState)),sigcorr_dg(in,:));
        case 2
            in = sum(ntrial_dg(:,2:4)>30,2)==3 & ~out;
            data(in,:,iState) = cellfun(@(x) nanmean(x(:,iState)),noisecorr_dg(in,:));
        case 3
            in = sum(nbin_sp(:,2:4)>50/binsize,2)==3 & ~out;
            data(in,:,iState) = cellfun(@(x) nanmean(x(:,iState)),noisecorr_sp(in,:));
    end
    msquare = NaN(4,4,nS);
    m = squeeze(data(:,:,iState));
    for i = 1:10
        msquare(x(i),y(i),:) = m(:,i);
        msquare(y(i),x(i),:) = m(:,i);
    end
    axes('Position',axpt(3,1,iState-1,1,axpt(1,10,1,1:9)));
    imagesc(squeeze(nanmean(msquare,3)));
    colormap(c);
    axis xy
    set(gca,'CLim',[-0.15 0.15],'FontSize',7,'LineWidth',0.35,'XTick',1:4,...
        'XTickLabel',{'Nomod';'Cluster 1';'Cluster 2';'Cluster 3'},'YTick',1:4,'YTickLabel',...
        {'Nomod';'Cluster 1';'Cluster 2';'Cluster 3'},...
        'XTickLabelRotation',45,'Box','off','TickDir','out');
    title(stateList{iState});
    if iState>2
        set(gca,'YTickLabel',[]);
    end
end
% cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig2');
% print(fHandle,'-depsc','-painters',...
%     ['correlation_',typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'.ai'])

%%
in = sum(nbin_sp(:,2:4)>50/binsize,2)==3 & ~out;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.15 20*0.18]);
hold on;
mwithin = squeeze(nanmean(data(in,x==y & x>1,2:4),2));
mbetween = squeeze(nanmean(data(in,x~=y & x>1 & y>1,2:4),2));
plot(1:3,mbetween,'Color',[0.8 0.8 0.8]);
plot(1:3,mwithin,'Color',[1 0.8 0.8]);
errorbar(1:3,nanmean(mbetween),nanstd(mbetween)/sqrt(size(mbetween,1)),'k')
errorbar(1:3,nanmean(mwithin),nanstd(mwithin)/sqrt(size(mwithin,1)),'r')
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',1:3,...
    'XTickLabel',{'low pupil';'high pupil';'mobile'},...
    'XTickLabelRotation',45,'YTick',-0.1:0.1:0.2)
ylim([-0.1 0.25])
xlim([0.5 3.5])
if iData==1
    ylabel('Signal correlation');
else
ylabel('Noise correlation')
end
cd('D:\OneDrive - UCSF\figures\2.allen\fig2');
print(fHandle,'-depsc','-painters',...
    ['correlation_',typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'_w_state.ai'])

tbl = simple_mixed_anova(cat(3,mwithin,mbetween));

%% for each cluster
titleList = {'iAct','dAct','Inh'};
[p,F] = deal(nan(3,3));
close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 7 20*0.175]);
for iClst = 1:3
    axes('Position',axpt(3,1,iClst,1,axpt(10,10,2:10,1:9)))
    hold on;
    mwithin = squeeze(nanmean(data(in,x==y & x==iClst+1,2:4),2));
    mbetween = squeeze(nanmean(data(in,x~=y & [(x==iClst+1 & y>1) | (x>1 & y==iClst+1)],2:4),2));
    plot(1:3,mbetween,'Color',[0.8 0.8 0.8]);
    plot(1:3,mwithin,'Color',[1 0.8 0.8]);
    errorbar(1:3,nanmean(mbetween),nanstd(mbetween)/sqrt(size(mbetween,1)),'k')
    errorbar(1:3,nanmean(mwithin),nanstd(mwithin)/sqrt(size(mwithin,1)),'r')
    set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',1:3,...
        'XTickLabel',{'Small pupil';'Large pupil';'Mobile'},...
        'XTickLabelRotation',45,'YTick',-0.15:0.15:0.3)
    ylim([-0.15 0.3])
    xlim([0.5 3.5])
    tbl = simple_mixed_anova(cat(3,mwithin,mbetween));
    text(0.7,-0.05,['p_c_l_u_s_t_e_r = ',num2str(tbl.pValue(5))],'FontSize',6);
    text(0.7,-0.09,['p_s_t_a_t_e = ',num2str(tbl.pValue(3))],'FontSize',6);
    text(0.7,-0.13,['p_c_x_s = ',num2str(tbl.pValue(7))],'FontSize',6);
    p(iClst,:) = tbl.pValue([5,3,7]);
    F(iClst,:) = tbl.F([5,3,7]);
    title(titleList{iClst});
    if iData==1
        ylabel('Signal correlation');
    else
        if iClst==1
        ylabel('Spontaneous correlation')
        else
            set(gca,'YTickLabel',[])
        end
    end
end
cd('D:\OneDrive - UCSF\figures\2.allen\fig2');
print(fHandle,'-depsc','-painters',...
    ['correlation_',typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'_w_state_eachcluster.ai'])

%%
data = NaN(nS,10);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.19 20/6]);
iState = 1;
switch iData
    case 1
        in = ~out;
        data(in,:) = cellfun(@(x) nanmean(x(:,iState)),sigcorr_dg(in,:));
    case 2
        in = ~out;
        data(in,:) = cellfun(@(x) nanmean(x(:,iState)),noisecorr_dg(in,:));
    case 3
        in = ~out;
        data(in,:) = cellfun(@(x) nanmean(x(:,iState)),noisecorr_sp(in,:));
end
msquare = NaN(4,4,nS);
for i = 1:10
    msquare(x(i),y(i),:) = data(:,i);
    msquare(y(i),x(i),:) = data(:,i);
end
axes('Position',axpt(10,10,3:10,1:9));
imagesc(squeeze(nanmean(msquare,3)));
colormap(c);
axis xy
set(gca,'CLim',[-0.15 0.15],'FontSize',7,'LineWidth',0.35,'XTick',1:4,...
    'XTickLabel',{'Nomod';'Cluster 1';'Cluster 2';'Cluster 3'},'YTick',1:4,'YTickLabel',...
    {'Nomod';'Cluster 1';'Cluster 2';'Cluster 3'},...
    'XTickLabelRotation',45,'Box','off','TickDir','out');

cd('D:\OneDrive - UCSF\figures\2.allen\fig2');
% print(fHandle,'-depsc','-painters',...
%     ['correlation_',typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'_all.ai'])

withinave = mean(data(:,x==y&x>1),2);
betweenave = mean(data(:,x~=y&x>1&y>1),2);
nomodave = data(:,x==y&x==1);
[~,p,~,stat] = ttest(withinave,betweenave);


%% across-ensemble correlation
close all
pairlist = [1,6,7,9];
in = ~out;
data = cellfun(@(x) mean(x(:,iState)),noisecorr_sp(in,:));
[p,t] = deal(nan(3,1));
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 20*0.22]);
plot(1:3,data(:,pairlist(2:end)),'Color',[0.6 0.6 0.6],'LineWidth',0.35);
hold on;
errorbar(1:3,mean(data(:,pairlist(2:end))),std(data(:,pairlist(2:end)))/sqrt(sum(in)),'Color','k','LineWidth',0.5);
plot([0 4],repmat(mean(data(:,pairlist(1))),1,2),'r--','LineWidth',0.35)
plot([0 4],[0 0],'k:')
for i = 1:3
[~,ptemp,~,stat] = ttest(data(:,pairlist(i+1)),data(:,pairlist(1)));
p(i) = ptemp*3;
t(i) = stat.tstat; 
if p(i)<0.05
text(i,0.13,'*')
end
end
xlim([0.5 3.5])
ylim([-0.2 0.2])
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',...
    1:3,'XTickLabel',{'iAct vs. dAct','iAct vs. Inh','dAct vs. Inh'},'XTickLabelRotation',45,...
    'YTick',-0.2:0.1:0.2)
ylabel('Spont. correlation')
cd('D:\OneDrive - UCSF\figures\2.allen\fig2');
print(fHandle,'-depsc','-painters',['correlation_',...
    typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'_acrossclusters.ai'])

%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 7 20*0.22]);
[p,F] = deal(nan(3,3));
for iClst = 1:3
    withincluster = x==y & x==iClst+1;
betweencluster = x~=y & [(x==iClst+1 & y>1) | (x>1 & y==iClst+1)];
    [c_withincl,c_betweencl] = deal(NaN(nS,2));
    for iS = find(in)'
        c_withincl(iS,1) = mean(cell2mat(cellfun(@(x,y) x(y,iState),noisecorr_sp(iS,withincluster)',...
            withinarea(iS,withincluster)','UniformOutput',false)));
        c_betweencl(iS,1) = mean(cell2mat(cellfun(@(x,y) x(y,iState),noisecorr_sp(iS,betweencluster)',...
            withinarea(iS,betweencluster)','UniformOutput',false)));
        c_withincl(iS,2) = mean(cell2mat(cellfun(@(x,y) x(~y,iState),noisecorr_sp(iS,withincluster)',...
            withinarea(iS,withincluster)','UniformOutput',false)));
        c_betweencl(iS,2) = mean(cell2mat(cellfun(@(x,y) x(~y,iState),noisecorr_sp(iS,betweencluster)',...
            withinarea(iS,betweencluster)','UniformOutput',false)));
    end
    axes('Position',axpt(3,1,iClst,1,axpt(10,10,2:10,1:9)))
    plot(1:4,[c_withincl,c_betweencl],'Color',[0.8 0.8 0.8]);
    hold on;
    errorbar(1:2,nanmean(c_withincl),nanstd(c_withincl)/sqrt(sum(in)),'r')
    errorbar(3:4,nanmean(c_betweencl),nanstd(c_betweencl)/sqrt(sum(in)),'k')
    xlim([0.5 4.5])
    ylim([-0.15 0.3])
    set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',1:4,...
        'XTickLabel',{'w/i area';'b/w areas';'w/i area';'b/w areas'},'XTickLabelRotation',45,...
        'YTick',-0.15:0.15:0.3);
    title(titleList{iClst})
    tbl = simple_mixed_anova(cat(3,c_withincl,c_betweencl));
    text(0.7,-0.02,['p_c_l_u_s_t_e_r = ',num2str(tbl.pValue(5))],'FontSize',6);
    text(0.7,-0.05,['p_a_r_e_a = ',num2str(tbl.pValue(3))],'FontSize',6);
    text(0.7,-0.08,['p_c_x_a = ',num2str(tbl.pValue(7))],'FontSize',6);
    p(iClst,:) = tbl.pValue([5,3,7]);
    F(iClst,:) = tbl.F([5,3,7]);
    if iClst==1
    ylabel('Noise correlation')
    else
        set(gca,'YTickLabel',[]);
    end
end
cd('D:\OneDrive - UCSF\figures\2.allen\fig2');
print(fHandle,'-depsc','-painters',['correlation_',...
    typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'_wbarea_eachcluster.ai'])
%% for each area
areaList = {'VISp';'VISl';'VISrl';'VISal';'VISpm';'VISam'};
nArea = length(areaList);
xx = [1,2,3,1,2,3];
yy = [1,1,1,2,2,2];
in = ~out;

%%
[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
structure = tag.info.structure(idx);
session = tag.info.session_id(idx);
sessionList = unique(session);
for iClst = 1:3
    for iA = 1:nArea
        for iS = 1:nS
        n(iA,iClst,iS) = sum(strcmp(structure,areaList{iA})&cluster_idx.vis{2}==iClst & session==sessionList(iS));        
        end
    end
end
for iA = 1:nArea
    subplot(2,3,iA)
    hold on;
    plot(1:3,squeeze(n(iA,:,:)),'Color',[0.8 0.8 0.8])
    errorbar(1:3,squeeze(mean(n(iA,:,:),3)),squeeze(std(n(iA,:,:),[],3))/sqrt(nS),'k')
    ylim([0 50]);
    xlim([0 4])
end
%%
invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis & celltype(:,iType)));
    unitid_tmp = T.unit_id(invis);
    [~,idx] = ismember(unitid_tmp,tag.info.unit_id);
    structure = tag.info.structure(idx);
    ccf = tag.info.ccf(idx,:);

    [in,idx] = ismember(unitid_tmp,unit_id.vis);
    clusteridx{iS} = zeros(length(idx),1);
    clusteridx{iS}(in) = cluster_idx.vis{k-1}(idx(in));





close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 7 20*0.22]);
[p,F] = deal(nan(3,3,nArea));
for iClst = 1:3
    withincluster = x==y & x==iClst+1;
    betweencluster = x~=y & [(x==iClst+1 & y>1) | (x>1 & y==iClst+1)];
    [c_withincl,c_betweencl,n_withincl,n_betweencl] = deal(NaN(nS,2,nArea));
    for iA = 1:nArea
        for iS = find(in)'
            c_withincl(iS,1,iA) = mean(cell2mat(cellfun(@(x,y,z) x(y&z(:,1)==areaList{iA},iState),...
                noisecorr_sp(iS,withincluster)',withinarea(iS,withincluster)',areas(iS,withincluster)', ...
                'UniformOutput',false)));
            c_betweencl(iS,1,iA) = mean(cell2mat(cellfun(@(x,y,z) x(y&z(:,1)==areaList{iA},iState), ...
                noisecorr_sp(iS,betweencluster)',withinarea(iS,betweencluster)',areas(iS,betweencluster)', ...
                'UniformOutput',false)));
            c_withincl(iS,2,iA) = mean(cell2mat(cellfun(@(x,y,z) x(~y&sum(z==areaList{iA},2)==1,iState), ...
                noisecorr_sp(iS,withincluster)',withinarea(iS,withincluster)',areas(iS,withincluster)', ...
                'UniformOutput',false)));
            c_betweencl(iS,2,iA) = mean(cell2mat(cellfun(@(x,y,z) x(~y&sum(z==areaList{iA},2)==1,iState), ...
                noisecorr_sp(iS,betweencluster)',withinarea(iS,betweencluster)',areas(iS,betweencluster)', ...
                'UniformOutput',false)));

            n_withincl(iS,:,iA) = sum(cell2mat(cellfun(@(x,y) [sum(x&sum(y==areaList{iA},2)>0),sum(~x&sum(y==areaList{iA},2)>0)],...
                withinarea(iS,withincluster)',areas(iS,withincluster)','UniformOutput',false)),1);
            n_betweencl(iS,:,iA) = sum(cell2mat(cellfun(@(x,y) [sum(x&sum(y==areaList{iA},2)>0),sum(~x&sum(y==areaList{iA},2)>0)],...
                withinarea(iS,betweencluster)',areas(iS,betweencluster)','UniformOutput',false)),1);
        end
        inanal = sum([squeeze(n_withincl(:,:,iA)),squeeze(n_betweencl(:,:,iA))]>3,2)==4;
        c_withincl(~inanal,:,iA) = nan;
        c_betweencl(~inanal,:,iA) = nan;

        axes('Position',axpt(3,1,iClst,1,axpt(3,2,xx(iA),yy(iA))))
        plot(1:4,[squeeze(c_withincl(:,:,iA)),squeeze(c_betweencl(:,:,iA))],'Color',[0.8 0.8 0.8]);
        hold on;
        errorbar(1:2,nanmean(squeeze(c_withincl(:,:,iA))),nanstd(squeeze(c_withincl(:,:,iA)))/sqrt(sum(inanal)),'r')
        errorbar(3:4,nanmean(squeeze(c_betweencl(:,:,iA))),nanstd(squeeze(c_betweencl(:,:,iA)))/sqrt(sum(inanal)),'k')
        xlim([0.5 4.5])
        ylim([-0.15 0.3])
        set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',1:4,...
            'XTickLabel',{'w/i area';'b/w areas';'w/i area';'b/w areas'},'XTickLabelRotation',45,...
            'YTick',-0.15:0.15:0.3);
        if yy(iA)==1
                title(titleList{iClst})
        end
                tbl = simple_mixed_anova(cat(3,squeeze(c_withincl(:,:,iA)),squeeze(c_betweencl(:,:,iA))));
%                 text(0.7,-0.02,['p_c_l_u_s_t_e_r = ',num2str(tbl.pValue(5))],'FontSize',6);
%                 text(0.7,-0.05,['p_a_r_e_a = ',num2str(tbl.pValue(3))],'FontSize',6);
%                 text(0.7,-0.08,['p_c_x_a = ',num2str(tbl.pValue(7))],'FontSize',6);
                p(iClst,:,iA) = tbl.pValue([5,3,7]);
                F(iClst,:,iA) = tbl.F([5,3,7]);
                if iClst==1 & xx(iA)==1
                    ylabel('Noise correlation')
                else
                    set(gca,'YTickLabel',[]);
                end
                if yy(iA)==1
                    set(gca,'XTickLabel',[]);
                end
                if iClst==1
                    text(1,0.28,areaList{iA},'FontSize',8)
                end

    end
end





%%
data = cellfun(@(x) x(:,iState),noisecorr_sp,'UniformOutput',false);
in = ~out;

withincluster = x==y & x>1;
betweencluster = x~=y & x>1 & y>1;
xwi = [0:100:500,1000];
ywi = NaN(nS,length(xwi)-1);

[c_withincl,c_betweencl] = deal(NaN(nS,2));
for iS = find(in)'
    c_withincl(iS,1) = mean(cell2mat(cellfun(@(x,y) x(y,iState),noisecorr_sp(iS,withincluster)',...
        withinarea(iS,withincluster)','UniformOutput',false)));
    c_betweencl(iS,1) = mean(cell2mat(cellfun(@(x,y) x(y,iState),noisecorr_sp(iS,betweencluster)',...
        withinarea(iS,betweencluster)','UniformOutput',false)));
    c_withincl(iS,2) = mean(cell2mat(cellfun(@(x,y) x(~y,iState),noisecorr_sp(iS,withincluster)',...
        withinarea(iS,withincluster)','UniformOutput',false)));
    c_betweencl(iS,2) = mean(cell2mat(cellfun(@(x,y) x(~y,iState),noisecorr_sp(iS,betweencluster)',...
        withinarea(iS,betweencluster)','UniformOutput',false)));
end

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.15 20*0.175]);
plot(1:4,[c_withincl,c_betweencl],'Color',[0.8 0.8 0.8]);
hold on;
errorbar(1:2,nanmean(c_withincl),nanstd(c_withincl)/sqrt(sum(in)),'r')
errorbar(3:4,nanmean(c_betweencl),nanstd(c_betweencl)/sqrt(sum(in)),'k')
xlim([0.5 4.5])
ylim([-0.1 0.2])
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',1:4,...
    'XTickLabel',{'w/i area';'b/w areas';'w/i area';'b/w areas'},'XTickLabelRotation',45,...
    'YTick',-0.1:0.1:0.3);
ylabel('Noise correlation')
cd('D:\OneDrive - UCSF\figures\2.allen\fig2');
print(fHandle,'-depsc','-painters',['correlation_',...
    typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'_wbarea.ai'])

tbl = simple_mixed_anova(cat(3,c_withincl,c_betweencl));


%% for each ensemble
data = cellfun(@(x) x(:,iState),noisecorr_sp,'UniformOutput',false);
in = ~out;
[p,F] = deal(nan(3,3));
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 7 20*0.175]);
for iClst = 1:3
    inpair = x==y & x==iClst+1;
    betweencluster = x~=y & [(x==iClst+1 & y>1) | (x>1 & y==iClst+1)];
    [c_withincl,c_betweencl] = deal(NaN(nS,2));
    for iS = find(in)'
        c_withincl(iS,1) = mean(cell2mat(cellfun(@(x,y) x(y,iState),noisecorr_sp(iS,inpair)',...
            withinarea(iS,inpair)','UniformOutput',false)));
        c_betweencl(iS,1) = mean(cell2mat(cellfun(@(x,y) x(y,iState),noisecorr_sp(iS,betweencluster)',...
            withinarea(iS,betweencluster)','UniformOutput',false)));
        c_withincl(iS,2) = mean(cell2mat(cellfun(@(x,y) x(~y,iState),noisecorr_sp(iS,inpair)',...
            withinarea(iS,inpair)','UniformOutput',false)));
        c_betweencl(iS,2) = mean(cell2mat(cellfun(@(x,y) x(~y,iState),noisecorr_sp(iS,betweencluster)',...
            withinarea(iS,betweencluster)','UniformOutput',false)));
    end
    axes('Position',axpt(3,1,iClst,1,axpt(10,10,2:10,1:9)))
    plot(1:4,[c_withincl,c_betweencl],'Color',[0.8 0.8 0.8]);
    hold on;
    errorbar(1:2,nanmean(c_withincl),nanstd(c_withincl)/sqrt(sum(in)),'r')
    errorbar(3:4,nanmean(c_betweencl),nanstd(c_betweencl)/sqrt(sum(in)),'k')
    xlim([0.5 4.5])
    ylim([-0.15 0.3])
    set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',1:4,...
        'XTickLabel',{'w/i area';'b/w areas';'w/i area';'b/w areas'},'XTickLabelRotation',45,...
        'YTick',-0.15:0.15:0.3);
    title(titleList{iClst})
    tbl = simple_mixed_anova(cat(3,c_withincl,c_betweencl));
    text(0.7,-0.02,['p_c_l_u_s_t_e_r = ',num2str(tbl.pValue(5))],'FontSize',6);
    text(0.7,-0.05,['p_a_r_e_a = ',num2str(tbl.pValue(3))],'FontSize',6);
    text(0.7,-0.08,['p_c_x_a = ',num2str(tbl.pValue(7))],'FontSize',6);
    p(iClst,:) = tbl.pValue([5,3,7]);
    F(iClst,:) = tbl.F([5,3,7]);
    if iClst==1
    ylabel('Noise correlation')
    else
        set(gca,'YTickLabel',[]);
    end
end
cd('D:\OneDrive - UCSF\figures\2.allen\fig2');
print(fHandle,'-depsc','-painters',['correlation_',...
    typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'_wbarea_eachcluster.ai'])

%%
iState = 1;
data = cellfun(@(x) x(:,iState),noisecorr_sp,'UniformOutput',false);
iClst = 2;
withincluster = x==y & x==iClst+1;
betweencluster = x~=y & [(x==iClst+1 & y>1) | (x>1 & y==iClst+1)];

% areas = cellfun(@(x) categorical(x),areas,'UniformOutput',false);
in = ~out;
[c_withincl,c_betweencl,n_withincl,n_betweencl] = deal(NaN(nS,nArea));
for iA = 1
    for iS = find(in)'
        c_withincl(iS,iA) = mean(cell2mat(cellfun(@(x,y,z) x(y&z(:,1)==areaList{iA},iState),...
            noisecorr_sp(iS,withincluster)',withinarea(iS,withincluster)',areas(iS,withincluster)', ...
            'UniformOutput',false)));
        c_betweencl(iS,iA) = mean(cell2mat(cellfun(@(x,y,z) x(y&z(:,1)==areaList{iA},iState), ...
            noisecorr_sp(iS,betweencluster)',withinarea(iS,betweencluster)',areas(iS,betweencluster)', ...
            'UniformOutput',false)));

        n_withincl(iS,iA) = length(cell2mat(cellfun(@(x,y,z) x(y&z(:,1)==areaList{iA},iState),...
            noisecorr_sp(iS,withincluster)',withinarea(iS,withincluster)',areas(iS,withincluster)', ...
            'UniformOutput',false)));
        n_betweencl(iS,iA) = length(cell2mat(cellfun(@(x,y,z) x(y&z(:,1)==areaList{iA},iState), ...
            noisecorr_sp(iS,betweencluster)',withinarea(iS,betweencluster)',areas(iS,betweencluster)', ...
            'UniformOutput',false)));
    end
end

%%
[p,tstat,nanal] = deal(nan(nArea,1));
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 9 20*0.19]);
for iA = 1:nArea
    axes('Position',axpt(nArea,1,iA,1,axpt(10,10,2:10,1:8),[0.025 0.05]))
    hold on;
    in = n_withincl(:,iA)>3 & n_betweencl(:,iA)>3;
    plot(1:2,[c_withincl(in,iA),c_betweencl(in,iA)],'Color',[0.8 0.8 0.8]);
    errorbar(1:2,nanmean([c_withincl(in,iA),c_betweencl(in,iA)]),...
        nanstd([c_withincl(in,iA),c_betweencl(in,iA)])/sqrt(sum(in)),'k')
    [~,p(iA),~,stat] = ttest(c_withincl(in,iA),c_betweencl(in,iA));
    tstat(iA) = stat.tstat;
    nanal(iA) = sum(in);
    if p(iA)<0.05
        text(1.5,0.18,'*')
    end
    ylim([-0.1 0.2])
    xlim([0.5 2.5])
    title(areaList{iA})
    set(gca,'XTick',[1:2],'XTickLabel',{'w/i cluster';'b/w clusters'},'XTickLabelRotation',45,...
        'FontSize',8,'LineWidth',0.35,'YTick',-0.1:0.1:0.2,'Box','off','TickDir','out');
    if iA>1
        set(gca,'YTickLabel',[])
    else
        ylabel('Spont. correlation')
    end
end
cd('D:\OneDrive - UCSF\figures\2.allen\fig2');
print(fHandle,'-depsc','-painters',['correlation_',...
    typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'_wiarea_dAct_eacharea.ai'])
