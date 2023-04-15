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
[clusteridx,distance,withinarea] = deal(cell(nS,1));

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
    load([sdir(sessionList(iS)),'_Events.mat'],'drifting_gratings','running_speed','pupil_data');
    load([sdir(sessionList(iS)),'_stimulus_evoked.mat'],'fr_gratings','beh_gratings');
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_ripples.mat'],'filtered_lfp','spontaneous_anal_win','spontaneous_win');

    % detect ripples
    win = [drifting_gratings.window(1,1) drifting_gratings.window(end,end)];
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
    [~,idx] = ismember(T.unit_id(invis),fr_gratings.unit_id);
    dglist = sort(unique(drifting_gratings.stimulus_condition_id));
    ndg = length(dglist);
    dgIndex = false(length(drifting_gratings.stimulus_condition_id),ndg);
    for idg = 1:ndg
        dgIndex(:,idg) = drifting_gratings.stimulus_condition_id==dglist(idg);
    end

    nCell = sum(invis);
    pupil_size = pupil_data.pupil_height/2.*pupil_data.pupil_width/2*pi;

    [sig_corr_dg,noise_corr_dg,noise_corr_sp] = deal(NaN(nCell*(nCell-1)/2,4));
    ripples = cell2mat(ripples);
    dg_ripple = cellfun(@(x) sum(ripples(:,1)>=x(1) & ripples(:,1)<=x(2)),...
        mat2cell(drifting_gratings.window,ones(size(drifting_gratings.window,1),1),2));
    for iState = 1:4
        [pairid,sig_corr_dg(:,iState),noise_corr_dg(:,iState),...
            ntrial_dg(iS,iState),pupil_dg(iS,iState)] =...
            paircorr(fr_gratings.spikeTime(idx),[0 2],binsize,dgIndex,...
            cell2mat(beh_gratings.velocity_ave'),speedlimit,unitid_tmp,iState,dg_ripple,...
            cell2mat(beh_gratings.pupil_area_ave'));

        [~,noise_corr_sp(:,iState),time_sp(iS,iState),pupil_sp(iS,iState),nbin_sp(iS,iState)] =...
            noisecorr_nostim(T.spike_time(idx),spontaneous_win,binsize,unitid_tmp,ripples(:,1),...
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

            sigcorr_dg{iS,kClst} = sig_corr_dg(in,:);
            noisecorr_dg{iS,kClst} = noise_corr_dg(in,:);
            noisecorr_sp{iS,kClst} = noise_corr_sp(in,:);

            distance{iS,kClst} = pair_distance(in);
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

stateList = {'all';'Immobile-low pupil';'Immobile-high pupil';'Mobile'};
dataList = {'dg_signal';'dg_noise';'sp_noise'};
c = cbrewer('div','RdBu',50);
c([1:5, 46:50],:) = [];

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
cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig2');
print(fHandle,'-depsc','-painters',...
    ['correlation_',typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'.ai'])

%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.19 20/6]);
iState = 1;
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
    axes('Position',axpt(10,10,3:10,1:9));
    imagesc(squeeze(nanmean(msquare,3)));
    colormap(c);
    axis xy
    set(gca,'CLim',[-0.15 0.15],'FontSize',7,'LineWidth',0.35,'XTick',1:4,...
        'XTickLabel',{'Nomod';'Cluster 1';'Cluster 2';'Cluster 3'},'YTick',1:4,'YTickLabel',...
        {'Nomod';'Cluster 1';'Cluster 2';'Cluster 3'},...
        'XTickLabelRotation',45,'Box','off','TickDir','out');

cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig2');
print(fHandle,'-depsc','-painters',...
    ['correlation_',typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'_all.ai'])


%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.15 20/6]);
hold on;
mwithin = squeeze(nanmean(data(:,x==y & x>1,2:4),2));
mbetween = squeeze(nanmean(data(:,x~=y & x>1 & y>1,2:4),2));
plot(1:3,mbetween,'Color',[0.8 0.8 0.8]);
plot(1:3,mwithin,'Color',[1 0.8 0.8]);
errorbar(1:3,nanmean(mbetween),nanstd(mbetween)/sqrt(size(mbetween,1)),'k')
errorbar(1:3,nanmean(mwithin),nanstd(mwithin)/sqrt(size(mwithin,1)),'r')
set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,'XTick',1:3,...
    'XTickLabel',{'low pupil';'high pupil';'mobile'},...
    'XTickLabelRotation',45,'YTick',-0.1:0.1:0.2)
ylim([-0.1 0.25])
xlim([0.5 3.5])
ylabel('Correlation')
print(fHandle,'-depsc','-painters',...
    ['correlation_',typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'_w_state.ai'])
%%

if iData==3
    iState =1;
    fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2 20*0.18]);
    hold on;
    
    mwithin = mean(squeeze(data(in,x==y & x>1,iState)),2);
    mbetween = mean(squeeze(data(in,x~=y & x>1 & y>1,iState)),2);
    
    plot(1:2,[mwithin mbetween],'Color',[0.8 0.8 0.8]);
    errorbar(1:2,nanmean([mwithin mbetween]),nanstd([mwithin,mbetween])/sqrt(sum(in)),'k','CapSize',3)
    [~,p(iData)] = ttest(mwithin,mbetween);
    set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',1:2,...
        'XTickLabel',{'w/i cluster';'b/w cluster'},'XTickLabelRotation',45,'YTick',0:0.1:0.2)
    ylim([-0.05 0.2])
    xlim([0.5 2.5])
    ylabel('Noise correlation')
    print(fHandle,'-depsc','-painters',...
        ['correlation_',typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'_wb.ai'])
    %%
    data = cellfun(@(x) x(:,iState),noisecorr_sp,'UniformOutput',false);
    
    withincluster = x==y & x>1;
    betweencluster = x~=y & x>1 & y>1;
    xwi = [0:100:500,1000];
    ywi = NaN(nS,length(xwi)-1);
    
    dwi = cellfun(@(x,y) x(y),distance,withinarea,'UniformOutput',false);
    cwi = cellfun(@(x,y) x(y),data,withinarea,'UniformOutput',false);
    [pvalwi,rwi] = deal(NaN(nS,2));
    betawi = NaN(nS,2,2);
    for iS = find(in)'
        [beta,~,stat] = glmfit(cell2mat(dwi(iS,betweencluster)')/10^3,cell2mat(cwi(iS,betweencluster)'));
        [~,~,b] = histcounts(cell2mat(dwi(iS,betweencluster)'),xwi);
        tempc = cell2mat(cwi(iS,betweencluster)');
        ywi(iS,:,2) = cellfun(@(x) nanmean(tempc(b==x)),num2cell(1:length(xwi)-1));
        betawi(iS,:,2) = beta;
        pvalwi(iS,2) = stat.p(2);
        rwi(iS,2) = corr(cell2mat(dwi(iS,betweencluster)'),cell2mat(cwi(iS,betweencluster)'));
        
        [beta,~,stat] = glmfit(cell2mat(dwi(iS,withincluster)')/10^3,cell2mat(cwi(iS,withincluster)'));
        [~,~,b] = histcounts(cell2mat(dwi(iS,withincluster)'),xwi);
        tempc = cell2mat(cwi(iS,withincluster)');
        ywi(iS,:,1) = cellfun(@(x) nanmean(tempc(b==x)),num2cell(1:length(xwi)-1));
        betawi(iS,:,1) = beta;
        pvalwi(iS,1) = stat.p(2);
        rwi(iS,1) = corr(cell2mat(dwi(iS,withincluster)'),cell2mat(cwi(iS,withincluster)'));
    end
end


iS = 3;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 2.5]);
hold on;
scatter(cell2mat(dwi(iS,betweencluster)'),cell2mat(cwi(iS,betweencluster)'),1,[0.6 0.6 0.6],'filled');
scatter(cell2mat(dwi(iS,withincluster)'),cell2mat(cwi(iS,withincluster)'),1,[1 0.6 0.6],'filled');
for i = 1:2
    if i==1
        incluster= betweencluster;
        clr = 'k';
    else
        incluster = withincluster;
        clr = 'r';
    end
    beta = glmfit(cell2mat(dwi(iS,incluster)'),cell2mat(cwi(iS,incluster)'));
    xx = [min(cell2mat(dwi(iS,incluster)')),max(cell2mat(dwi(iS,incluster)'))];
    plot(xx,beta(1)+beta(2)*xx,'Color',clr);
    text(0,0.8-(i-1)*0.12,['y = ',num2str(round(betawi(iS,2,i)*100)/100),...
        'x+',num2str(round(betawi(iS,1,i)*100)/100)],'FontSize',4,'Color',clr);
end
set(gca,'YLim',[-0.6 0.9],'XTick',0:200:600,'XTickLabel',0:0.2:0.6,...
    'YTick',-0.6:0.3:0.9,'Box','off','TickDir','out','FontSize',5,'XLim',[-50 600])
title('Within area')
ylabel('Correlation')
xlabel('Distance (mm)')
print(fHandle,'-dtiff','-r600',['correlation_distance_',...
    typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'_ex.tif'])

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 2.5]);
for i = 1:2
    subplot(1,2,2/i)
    plot(1:2,squeeze(betawi(:,i,:)),'Color',[0.6 0.6 0.6]);
    hold on;
    errorbar(1:2,nanmean(squeeze(betawi(:,i,:))),nanstd(squeeze(betawi(:,i,:)))/...
        sqrt(sum(~isnan(betawi(:,i,1)))),'k');
    plot([0.5 2.5],[0 0],'k:');
    xlim([0.5 2.5])
    if i==1
        ylim([-0.1 0.3])
        set(gca,'YTick',-0.1:0.1:0.3);
        ylabel('Y intercept');
    else
        ylim([-1 0.5])
        set(gca,'YTick',-1:0.5:0.5);
        ylabel('Slope');
    end
    set(gca,'XTick',[1,2],'XTickLabel',{'within';'between'},'Box','off','TickDir','out',...
        'FontSize',5,'LineWidth',0.35,'XTickLabelRotation',45);
end
print(fHandle,'-dtiff','-r600',['correlation_distance_',...
    typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'_fitting.tif'])

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

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 2.5]);
plot(1:4,[c_withincl,c_betweencl],'Color',[0.8 0.8 0.8]);
hold on;
errorbar(1:2,nanmean(c_withincl),nanstd(c_withincl)/sqrt(sum(in)),'r')
errorbar(3:4,nanmean(c_betweencl),nanstd(c_betweencl)/sqrt(sum(in)),'k')
xlim([0.5 4.5])
ylim([-0.1 0.2])
set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,'XTick',1:4,...
    'XTickLabel',{'w/in area';'b/w areas';'w/in area';'b/w areas'},'XTickLabelRotation',45,...
    'YTick',-0.1:0.1:0.3);
ylabel('Correlation')
print(fHandle,'-dtiff','-r600',['correlation_',...
    typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'_wbarea.tif'])

%%
function [pairid,signalcorr,noisecorr,ntrial,pupil] =...
    paircorr(spikeTime,win,binsize,condition,speed,speedlimit,unit_id,state,trialripple,pupilsize)


spkhist = cellfun(@(x) cell2mat(cellfun(@(y) histcounts(y,win(1):binsize:win(2)),x,...
    'UniformOutput',false)),spikeTime,'UniformOutput',false);

nCell = length(spikeTime);
out = ones(nCell,nCell);
out = triu(out);
out = out(:);

x = repmat(unit_id,1,length(unit_id));
y = repmat(unit_id',length(unit_id),1);

pairid = [x(:),y(:)];
pairid = pairid(~out,:);

nCondition = size(condition,2);
[signal,noise] = deal(cell(nCell,nCondition));

switch state
    case 1 % no ripple
        instate = trialripple==0;
    case {2,3}% immobile & no ripple & low (or high) pupil (bottom (or top) 50% pupil within immobile&no ripple condition)
        pupilmobile = nanmean(pupilsize(abs(speed)>speedlimit));
        inconditiontemp = abs(speed)<speedlimit & trialripple==0;
        if state==2
            instate = inconditiontemp & pupilsize<pupilmobile*0.5;
        elseif state==3
            instate = inconditiontemp & pupilsize>=pupilmobile*0.5;
        end
    case 4 % mobile
        instate = abs(speed)>speedlimit;
end
if sum(instate)>0

for iC = 1:nCondition
    incondition = instate & condition(:,iC);
    if sum(incondition)==0
        continue;
    end
    signal(:,iC) = cellfun(@(x) nanmean(x(incondition,:)),spkhist,'UniformOutput',false);
    noise(:,iC) = cellfun(@(x,y) x(incondition,:)-repmat(y,sum(incondition),1),spkhist,signal(:,iC),'UniformOutput',false);
    noise(:,iC) = cellfun(@(x) x(:)',noise(:,iC),'UniformOutput',false);
end
signal = cell2mat(signal);
scorr = corr(signal');
scorr = scorr(:);
signalcorr = scorr(~out);

noise = cell2mat(noise);
ncorr = corr(noise');
ncorr = ncorr(:);
noisecorr = ncorr(~out);
ntrial = sum(instate);
pupil = nanmean(pupilsize(instate));
else
    ntrial = 0;
    pupil = NaN;
    [signalcorr,noisecorr] = deal(NaN(size(pairid,1),1));
end

end

function [pairid,noisecorr,avepupil,avespeed,nbin] = noisecorr_nostim(spikeTime,win,binsize,unit_id,rippletime,...
    iState,pupiltime,pupilsize,speedtime,speed,speedlimit)

spiketmp = cellfun(@(x) x(x>=win(1) & x<=win(2)),spikeTime,'UniformOutput',false);
spkhist = cell2mat(cellfun(@(x) histcounts(x,win(1):binsize:win(2)),spiketmp,'UniformOutput',false));

rippletmp = rippletime(rippletime>=win(1) & rippletime<=win(2));
ripplehist = histcounts(rippletmp,win(1):binsize:win(2));

[~,~,bin] = histcounts(pupiltime,win(1):binsize:win(2));
pupilhist = cellfun(@(x) nanmean(pupilsize(bin==x)),num2cell(1:floor(diff(win)/binsize)));

[~,~,bin] = histcounts(speedtime,win(1):binsize:win(2));
speedhist = cellfun(@(x) nanmean(speed(bin==x)),num2cell(1:floor(diff(win)/binsize)));

switch iState
    case 1 % every state & no ripple
        incondition = find(ripplehist==0);
    case {2,3} % immobile & no ripple & low-high level of pupil
        pupilmobile = nanmean(pupilhist(abs(speedhist)>speedlimit));
        inconditiontemp = abs(speedhist)<speedlimit & ripplehist==0;
        if iState==2
            incondition = find(inconditiontemp & pupilhist<0.5*pupilmobile);
        elseif iState==3
            incondition = find(inconditiontemp & pupilhist>=0.5*pupilmobile);
        end        
    case 4 % mobile
        incondition = find(abs(speedhist)>speedlimit);
end

if ~isempty(incondition)
    spkhist = spkhist(:,incondition);
    
    nCell = length(spikeTime);
    out = ones(nCell,nCell);
    out = triu(out);
    out = out(:);
    
    x = repmat(unit_id,1,length(unit_id));
    y = repmat(unit_id',length(unit_id),1);
    
    pairid = [x(:),y(:)];
    pairid = pairid(~out,:);
    
    noisecorr = corr(spkhist');
    noisecorr = noisecorr(:);
    noisecorr = noisecorr(~out);
    
    nbin = length(incondition);
    
    avepupil = nanmean(pupilhist(incondition));
    avespeed = nanmean(speedhist(incondition));
else
    nbin = 0;
    [avepupil,avespeed] = deal(NaN);
    noisecorr = NaN(size(pairid,1),1);
end
end