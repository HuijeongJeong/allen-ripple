clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');

k = 3;

unitid = unit_id.vis;
[~,idx] = ismember(unitid,tag.info.unit_id);
sessionid = tag.info.session_id(idx);
sessionList = unique(sessionid);
nS = length(sessionList);
celltype = [tag.celltype.rs, tag.celltype.fs];

[sigcorr_dg,noisecorr_dg,noisecorr_sp] = deal(cell(nS,k*(k+1)/2+k+1));
[ntrial_dg,time_sp,pupil_dg,pupil_sp] = deal(NaN(nS,3,2));
[structure,clusteridx,ccf,distance] = deal(cell(nS,1));

speedlimit = 2;
binsize = 2;

threshold = [2 5];   %threshold for ripple beginning/end and peak
duration = [30 20 250];    %min inter-ripple inter, min ripple length, max ripple length
Fs = 1250;

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
    structure{iS} = tag.info.structure(idx);
    ccf{iS} = tag.info.ccf(idx,:);

    [in,idx] = ismember(unitid_tmp,unitid);
    clusteridx{iS} = zeros(length(idx),1);
    clusteridx{iS}(in) = cluster_idx.vis{k-1}(idx(in));

    [~,idx] = ismember(T.unit_id(invis),fr_gratings.unit_id);
    dglist = sort(unique(drifting_gratings.stimulus_condition_id));
    ndg = length(dglist);
    dgIndex = false(length(drifting_gratings.stimulus_condition_id),ndg);
    for idg = 1:ndg
        dgIndex(:,idg) = drifting_gratings.stimulus_condition_id==dglist(idg);
    end

    nCell = sum(invis);

    win_sp = cell(2,1);
    spontaneous_mobile = [spontaneous_win(1), spontaneous_anal_win(1,1);...
        spontaneous_anal_win(1:end-1,2), spontaneous_anal_win(2:end,1);...
        spontaneous_anal_win(end,1) spontaneous_win(2)];
    spontaneous_mobile(diff(spontaneous_mobile,[],2)<=binsize,:) = [];
    win_sp{1} = spontaneous_anal_win;
    win_sp{2} = spontaneous_anal_win;
    win_sp{3} = spontaneous_mobile;

    pupil_size = pupil_data.pupil_height/2.*pupil_data.pupil_width/2*pi;
    starts = find(diff([1;running_speed.immobile])==-1);
    stops = find(diff([running_speed.immobile;1])==1);
    mobile_win = running_speed.time([starts,stops]);
    inmobile = sum(cell2mat(cellfun(@(x) pupil_data.time>=x(1) & pupil_data.time<=x(2),...
        mat2cell(mobile_win,ones(size(mobile_win,1),1),2),'UniformOutput',false)'),2)>0 &...
        pupil_data.time>=drifting_gratings.window(1,1) & pupil_data.time<=drifting_gratings.window(end,end);

    mobilepupil = pupil_size(inmobile);
    mobilepupil = nanmedian(mobilepupil);
    normpupil = pupil_size/mobilepupil;

    [sig_corr_dg,noise_corr_dg,noise_corr_sp] = deal(NaN(nCell*(nCell-1)/2,3));

    ripples = cell2mat(ripples);
    dg_ripple = cellfun(@(x) sum(ripples(:,1)>=x(1) & ripples(:,1)<=x(2)),...
        mat2cell(drifting_gratings.window,ones(size(drifting_gratings.window,1),1),2));
    for iState = 1:3
        [pairid,sig_corr_dg(:,iState),noise_corr_dg(:,iState),...
            ntrial_dg(iS,iState),pupil_dg(iS,iState)] =...
            paircorr(fr_gratings.spikeTime(idx),[0 2],binsize,dgIndex,...
            cell2mat(beh_gratings.velocity_ave'),speedlimit,unitid_tmp,iState,dg_ripple,...
            cell2mat(beh_gratings.pupil_area_ave')/mobilepupil);

        [~,noise_corr_sp(:,iState),time_sp(iS,iState),pupil_sp(iS,iState)] =...
            noisecorr_nostim(T.spike_time(idx),win_sp{iState},binsize,unitid_tmp,ripples(:,1),...
            iState,pupil_data.time,normpupil);
    end
    ccfout = sum(ccf{iS}<0,2)>0;
    ccf{iS}(ccfout,:) = NaN;
    pair_distance = pdist(ccf{iS})';

    [~,idx] = ismember(pairid,unitid_tmp);
    cidxtmp = [clusteridx{iS}(idx(:,1)),clusteridx{iS}(idx(:,2))];
    kClst = 1;
    for iClst = 1:k+1
        for jClst = iClst:k+1
            in = cidxtmp(:,1)==iClst-1 & cidxtmp(:,2)==jClst-1;

            sigcorr_dg{iS,kClst} = sig_corr_dg(in,:);
            noisecorr_dg{iS,kClst} = noise_corr_dg(in,:);
            noisecorr_sp{iS,kClst} = noise_corr_sp(in,:);

            distance{iS,kClst} = pair_distance(in);
            kClst = kClst+1;
        end
    end
end

nmod = cell2mat(cellfun(@(x) [sum(x==1),sum(x==2),sum(x==3)],clusteridx,'UniformOutput',false));
out = sum(nmod==0,2)>0;

clr = {[0.6 0.6 0.6],[0 0 0]; [0.6 0.6 1],[0 0 1]; [1 0.6 0.6],[1 0 0]};

[panova_noise,panova_signal] = deal(NaN(4,2));
[c_noise,c_signal] = deal(NaN(4,2,45));

xticks = {'Nomod vs. Nomod';'Nomod vs. Iact';'Nomod vs. Dact';'Nomod vs. Inh';...
    'Iact vs. Iact';'Iact vs. Dact';'Iact vs. Inh';'Dact vs. Dact';'Dact vs. Inh';'Inh vs. Inh'};
refList = {'d';'i'};

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 6]);
axes('Position',axpt(2,2,1,1,axpt(1,10,1,1:9)));
hold on;
for iState = 1:3
    data = cellfun(@(x) nanmean(x(:,iState)),sigcorr_dg(~out,:));
    plot([1:4,6:8,10:11,13],data,'Color',clr{iState,1},'LineWidth',0.3);
end
data = cell(3,1);
for iState = 1:3
    data{iState} = cellfun(@(x) nanmean(x(:,iState)),sigcorr_dg(~out,:));
    errorbar([1:4,6:8,10:11,13],nanmean(data{iState}),nanstd(data{iState})/sqrt(size(data{iState},1)),...
        'Color',clr{iState,2},'LineWidth',0.5,'CapSize',3)
    [panova_signal(1,iState),~,stat] = anova1(data{iState},[],'off');
    ctmp = multcompare(stat,'display','off');
    c_signal(1,iState,:) = ctmp(:,6);
end
xlim([0 14])
ylim([-0.2 0.4]);
plot([0 14],[0 0],' k:','LineWidth',0.35);
set(gca,'XTick',[1:4,6:8,10:11,13],'XTickLabel',[],'Box','off','TickDir',...
    'out','FontSize',4,'LineWidth',0.35,'YTick',-0.2:0.2:0.4);
ylabel('Signal correlation','FontSize',5);
title('Drifting gratings','FontSize',5);

for iData = 1:2
    axes('Position',axpt(2,2,iData,2,axpt(1,10,1,1:9)));
    hold on;
    data = [];
    for iState = 1:3
        switch iData
            case 1
                data = cellfun(@(x) nanmean(x(:,iState)),noisecorr_dg(~out,:));
            case 2
                data = cellfun(@(x) nanmean(x(:,iState)),noisecorr_sp(~out,:));
        end
        plot([1:4,6:8,10:11,13],data,'Color',clr{iState,1},'Linewidth',0.3);
    end
    
    data = cell(3,1);
    for iState = 1:3
        switch iData
            case 1
                data{iState} = cellfun(@(x) nanmean(x(:,iState)),noisecorr_dg(~out,:));
            case 2
                data{iState} = cellfun(@(x) nanmean(x(:,iState)),noisecorr_sp(~out,:));
        end
        errorbar([1:4,6:8,10:11,13],nanmean(data{iState}),nanstd(data{iState})/sqrt(size(data{iState},1)),...
            'Color',clr{iState,2},'LineWidth',0.5,'CapSize',3)
        [panova_noise(iData,iState),~,stat] = anova1(data{iState},[],'off');
        ctmp = multcompare(stat,'display','off');
        c_noise(iData,iState,:) = ctmp(:,6);
    end
    
    xlim([0 14])
    ylim([-0.2 0.4]);
    plot([0 14],[0 0],' k:','LineWidth',0.35);
    set(gca,'XTick',[1:4,6:8,10:11,13],'XTickLabel',xticks,'XTickLabelRotation',45,...
        'YTick',-0.2:0.2:0.4,'Box','off','TickDir','out','FontSize',4,'LineWidth',0.35);
    if iData==1
        ylabel('Noise correlation','FontSize',5);
    else
        set(gca,'YTickLabel',[]);
        title('Spontaneous','FontSize',5);
    end
end
print(fHandle,'-dtiff','-r600',...
    ['D:\OneDrive\Research\DataFig\1.allen-andermann\clustering\corr\pair_correlation_',...
    typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple.tif']);

    
%%
msquare = cell(3,3);
x = [1 1 1 1 2 2 2 3 3 4];
y = [1 2 3 4 2 3 4 3 4 4];
for iData = 1:3
    for iState = 1:3
        switch iData
            case 1
                data = cellfun(@(x) nanmean(x(:,iState)),sigcorr_dg(~out,:));
            case 2
                data = cellfun(@(x) nanmean(x(:,iState)),noisecorr_dg(~out,:));
            case 3
                data = cellfun(@(x) nanmean(x(:,iState)),noisecorr_sp(~out,:));
        end
        m = nanmean(data);
        msquare{iData,iState} = NaN(4,4);
        for i = 1:10
            msquare{iData,iState}(x(i),y(i)) = m(i);
            msquare{iData,iState}(y(i),x(i)) = m(i);
        end
    end
end

%%
jS = [2 1 3];
stateList = {'Immobile w/o ripple';'Immobile with ripple';'Mobile'};
dataList = {'Gratings-signal';'Gratnig-noise';'Spontaneous'};
c = cbrewer('div','RdBu',50);
c([1:5, 46:50],:) = [];
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 8]);
for iData = 1:3
    for iState = 1:3
        axes('Position',axpt(3,3,jS(iState),iData,axpt(15,10,2:15,1:9)));
        imagesc(msquare{iData,iState});
        colormap(c);
        axis xy
        set(gca,'CLim',[-0.15 0.15],'FontSize',5,'LineWidth',0.35,'XTick',1:4,...
            'XTickLabel',{'Nomod';'I act';'D act';'Inh'},'YTick',1:4,'YTickLabel',...
            {'Nomod';'I act';'D act';'Inh'},'XTickLabelRotation',45,'Box','off','TickDir','out');
        if iData<3
                        set(gca,'XTickLabel',[]);

        if iData==1
            title(stateList{iState},'FontSize',5,'Color',clr{iState,2});
        end
        end
        if jS(iState)>1
            set(gca,'YTickLabel',[]);
        else
            ylabel(dataList{iData},'FontSize',5);
        end
    end
end
    print(fHandle,'-dtiff','-r600',...
        ['D:\OneDrive\Research\DataFig\1.allen-andermann\clustering\corr\correlation_',...
        typeList{iType},'_',num2str(binsize*1000),'bin_wo_ripple.tif'])
    
    
refList = {'Dorsal';'Intermediate'};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 7 3]);
for iRef = 1:2
    axes('Position',axpt(2,1,iRef,1,axpt(2,12,1,2:10,axpt(15,1,2:15,1),[0.2 0.05])))
    plot(1:3,ntrial_dg(~out,:,iRef),'Color',[0.8 0.8 0.8]);
    hold on;
    errorbar(1:3,nanmean(squeeze(ntrial_dg(~out,:,iRef))),nanstd(squeeze(ntrial_dg(~out,:,iRef)))/sqrt(sum(~out)),'k');
    xlim([0 4])
    ylim([0 500]);
    set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,'XTick',1:3,...
        'YTick',0:100:500,'XTickLabel',{'Imb w/o rp';'Imb w/ rp';'Mobile'},...
        'XTickLabelRotation',45);
    if iRef==2
        set(gca,'YTickLabel',[]);
    else
        ylabel({'Drifting gratings';'Trial'},'FontSize',5);
    end
    title(refList{iRef},'FontSize',5);

    axes('Position',axpt(2,1,iRef,1,axpt(2,12,2,2:10,axpt(15,1,2:15,1),[0.2 0.05])))
    plot(1:3,time_sp(~out,:,iRef),'Color',[0.8 0.8 0.8]);
    hold on;
    errorbar(1:3,nanmean(squeeze(time_sp(~out,:,iRef))),nanstd(squeeze(time_sp(~out,:,iRef)))/sqrt(sum(~out)),'k');
    xlim([0 4])
    ylim([0 1500]);
    set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,'XTick',1:3,...
        'YTick',0:500:1500,'XTickLabel',{'Imb w/o rp';'Imb w/ rp';'Mobile'},...
        'XTickLabelRotation',45);
    if iRef==2
        set(gca,'YTickLabel',[]);
    else
        ylabel({'Spontaneous';'Time (s)'},'FontSize',5);
    end
    title(refList{iRef},'FontSize',5);
end
print(fHandle,'-dtiff','-r600','D:\heejeong\OneDrive\Research\DataFig\1.allen-andermann\clustering\corr\analysis_state_compare.tif');

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 7 3]);
for iRef = 1:2
    axes('Position',axpt(2,1,iRef,1,axpt(2,12,1,2:10,axpt(15,1,2:15,1),[0.2 0.05])))
    plot(1:3,pupil_dg(~out,:,iRef),'Color',[0.8 0.8 0.8]);
    hold on;
    errorbar(1:3,nanmean(squeeze(pupil_dg(~out,:,iRef))),nanstd(squeeze(pupil_dg(~out,:,iRef)))/sqrt(sum(~out)),'k');
    xlim([0 4])
    ylim([0 1.5]);
    set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,'XTick',1:3,...
        'YTick',0:0.5:1.5,'XTickLabel',{'Imb w/o rp';'Imb w/ rp';'Mobile'},...
        'XTickLabelRotation',45);
    if iRef==2
        set(gca,'YTickLabel',[]);
    else
        ylabel({'Drifting gratings';'Normalized pupil size'},'FontSize',5);
    end
    title(refList{iRef},'FontSize',5);

    axes('Position',axpt(2,1,iRef,1,axpt(2,12,2,2:10,axpt(15,1,2:15,1),[0.2 0.05])))
    plot(1:3,pupil_sp(~out,:,iRef),'Color',[0.8 0.8 0.8]);
    hold on;
    errorbar(1:3,nanmean(squeeze(pupil_sp(~out,:,iRef))),nanstd(squeeze(pupil_sp(~out,:,iRef)))/sqrt(sum(~out)),'k');
    xlim([0 4])
    ylim([0 1.5]);
    set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,'XTick',1:3,...
        'YTick',0:0.5:1.5,'XTickLabel',{'Imb w/o rp';'Imb w/ rp';'Mobile'},...
        'XTickLabelRotation',45);
    if iRef==2
        set(gca,'YTickLabel',[]);
    else
        ylabel({'Spontaneous';'Normalized pupil size'},'FontSize',5);
    end
    title(refList{iRef},'FontSize',5);
end
print(fHandle,'-dtiff','-r600','D:\heejeong\OneDrive\Research\DataFig\1.allen-andermann\clustering\corr\analysis_state_pupil_compare.tif');

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
for iC = 1:nCondition
    switch state
        case 1
            incondition = condition(:,iC) & abs(speed)<speedlimit & trialripple==0;
        case 2
            incondition = condition(:,iC) & abs(speed)<speedlimit & trialripple>0;
        case 3
            incondition = condition(:,iC) & abs(speed)>speedlimit;
    end
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
switch state
    case 1
        ntrial = sum(abs(speed)<speedlimit & trialripple==0);
        pupil = nanmean(pupilsize(abs(speed)<speedlimit & trialripple==0));
    case 2
        ntrial = sum(abs(speed)<speedlimit & trialripple>0);
        pupil = nanmean(pupilsize(abs(speed)<speedlimit & trialripple>0));
    case 3
        ntrial = sum(abs(speed)>speedlimit);
        pupil = nanmean(pupilsize(abs(speed)>speedlimit));
end
end



function [pairid,noisecorr,time,pupil] = noisecorr_nostim(spikeTime,win,binsize,unit_id,rippletime,iRipple,pupiltime,pupilsize)

nWin = size(win,1);
[spkhist,ripplehist,pupilhist] = deal(cell(1,nWin));
for iW = 1:nWin
    if diff(win(iW,:))<binsize
        continue;
    end
    spiketmp = cellfun(@(x) x(x>=win(iW,1) & x<=win(iW,2)),spikeTime,'UniformOutput',false);
    spkhist{iW} = cell2mat(cellfun(@(x) histcounts(x,win(iW,1):binsize:win(iW,2)),spiketmp,'UniformOutput',false));

    rippletmp = rippletime(rippletime>=win(iW,1) & rippletime<=win(iW,2));
    ripplehist{iW} = histcounts(rippletmp,win(iW,1):binsize:win(iW,2));

    [~,~,bin] = histcounts(pupiltime,win(iW,1):binsize:win(iW,2));
    pupilhist{iW} = cellfun(@(x) nanmean(pupilsize(bin==x)),num2cell(1:floor(diff(win(iW,:))/binsize)));
end

spkhist = cell2mat(spkhist);
ripplehist = cell2mat(ripplehist);
pupilhist = cell2mat(pupilhist);
switch iRipple
    case 1
        spkhist(:,ripplehist>0) = [];
        pupilhist(ripplehist>0) = [];
    case 2
        spkhist(:,ripplehist==0) = [];
        pupilhist(ripplehist==0) = [];
end

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

time = length(spkhist)*binsize;
pupil = nanmean(pupilhist);
end