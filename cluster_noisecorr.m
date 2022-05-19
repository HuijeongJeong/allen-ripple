clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');

%% setting up 
k = 3; % number of cluster 

unitid = unit_id.vis;
[~,idx] = ismember(unitid,tag.info.unit_id);
sessionid = tag.info.session_id(idx);
sessionList = unique(sessionid);
nS = length(sessionList);
celltype = [tag.celltype.rs, tag.celltype.fs];

% noise correlation during spontaneous period
% (session, cluster pair, celltype) 
% cluster pair: nomod-nomod, nomod-Iact, nomod-Dact, nomod-Inh, Iact-Dact,...
% cell type: all pairs, rs vs. rs, fs vs. fs, rs vs. fs
[noisecorr_sp,sigcorr_dg] = deal(cell(nS,k*(k+1)/2+k+1,4)); 

% noise correlation during spontaneous period within each area or across areas
% (session, cluster pair, within/across areas, celltype)
[noisecorr_sp_area,sigcorr_dg_area] = deal(cell(nS,k*(k+1)/2+k+1,2,4)); 

threshold = [2 5];   %threshold for ripple beginning/end and peak
duration = [30 20 250];    %min inter-ripple inter, min ripple length, max ripple length
Fs = 1250;

binsize = 2; % unit: second
speedlimit = 2;
%% 
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_Events.mat'],'running_speed','pupil_data','drifting_gratings');
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_ripples.mat'],'filtered_lfp','spontaneous_anal_win','spontaneous_CA1_ripple');
    load([sdir(sessionList(iS)),'_stimulus_evoked.mat'],'fr_gratings','beh_gratings');
    
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
    
    
    invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis));
    unitid_tmp = T.unit_id(invis);
    [~,idx] = ismember(unitid_tmp,tag.info.unit_id);
    
    structure = categorical(tag.info.structure(idx)); % structure acronym
    
    [r,c] = find(celltype(idx,:));
    [r,sortidx] = sort(r);
    ctype = zeros(sum(invis),1);
    ctype(r) = c(sortidx); % cell type: 1 regular-spiking, 2 fast-spiking
    
    [in,idx] = ismember(unitid_tmp,unitid);
    clusteridx = zeros(length(idx),1);
    clusteridx(in) = cluster_idx.vis{k-1}(idx(in)); % cluster idx: 1 Iact, 2 Dact, 3 Inh
    
    % if total numer of ripple-modulated units are less than 5, exclude the session from analysis
    if sum(in)<5
        continue; 
    end
    
    [~,idx] = ismember(T.unit_id(invis),fr_gratings.unit_id);
    dglist = sort(unique(drifting_gratings.stimulus_condition_id));
    ndg = length(dglist);
    dgIndex = false(length(drifting_gratings.stimulus_condition_id),ndg);
    for idg = 1:ndg
        dgIndex(:,idg) = drifting_gratings.stimulus_condition_id==dglist(idg);
    end
    
    pupil_size = pupil_data.pupil_height/2.*pupil_data.pupil_width/2*pi;
    mobile_win = running_speed.time([find(diff([1;running_speed.immobile])==-1),...
        find(diff([running_speed.immobile;1])==1)]);
    inmobile = sum(cell2mat(cellfun(@(x) pupil_data.time>=x(1) & pupil_data.time<=x(2),...
        mat2cell(mobile_win,ones(size(mobile_win,1),1),2),'UniformOutput',false)'),2)>0 &...
        pupil_data.time>=spontaneous_anal_win(1,1) & pupil_data.time<=spontaneous_anal_win(end,2);
    mobilepupil = pupil_size(inmobile);
    mobilepupil = nanmedian(mobilepupil); % median of pupil size during mobile state
    normpupil = pupil_size/mobilepupil; % normalized pupil size by mobile pupil size
    
    rippletime = cell2mat(spontaneous_CA1_ripple.ripple(...
        ~cellfun(@isempty,spontaneous_CA1_ripple.relative_location)));
    rippletime = rippletime(:,1); % all ripples from dorsal & intermediate ca1
    
    nCell = sum(invis);
    [pairid,noise_corr_sp,time_sp,pupil_sp] =...
        noisecorr_nostim(T.spike_time(invis),spontaneous_anal_win,binsize,...
        unitid_tmp,rippletime,1,pupil_data.time,normpupil);
    
    ripples = cell2mat(ripples);
    dg_ripple = cellfun(@(x) sum(ripples(:,1)>=x(1) & ripples(:,1)<=x(2)),...
        mat2cell(drifting_gratings.window,ones(size(drifting_gratings.window,1),1),2));
    [~,sig_corr_dg,noise_corr_dg] =...
            paircorr(fr_gratings.spikeTime(idx),[0 2],binsize,dgIndex,...
            cell2mat(beh_gratings.velocity_ave'),speedlimit,unitid_tmp,1,dg_ripple,...
            cell2mat(beh_gratings.pupil_area_ave')/mobilepupil);

    [~,idx] = ismember(pairid,unitid_tmp);
    cidxtmp = [clusteridx(idx(:,1)),clusteridx(idx(:,2))];
    ctypetmp = [ctype(idx(:,1)),ctype(idx(:,2))];
    area = [structure(idx(:,1)), structure(idx(:,2))];
    
    kClst = 1;
    for iClst = 1:k+1
        for jClst = iClst:k+1
            in = cidxtmp(:,1)==iClst-1 & cidxtmp(:,2)==jClst-1;
            for i = 1:4
                if i>1 & i<4
                    inct = sum(ctypetmp==i-1,2)==2;
                elseif i==4
                    inct = sum(ctypetmp==1,2)==1 & sum(ctypetmp==2,2)==1;
                elseif i==1
                    inct = in;
                end
                noisecorr_sp{iS,kClst,i} = noise_corr_sp(in&inct);
                noisecorr_sp_area{iS,kClst,1,i} = noise_corr_sp(in & inct & area(:,1)==area(:,2));
                noisecorr_sp_area{iS,kClst,2,i} = noise_corr_sp(in & inct & area(:,1)~=area(:,2));
                
                sigcorr_dg{iS,kClst,i} = sig_corr_dg(in&inct);
                sigcorr_dg_area{iS,kClst,1,i} = sig_corr_dg(in & inct & area(:,1)==area(:,2));
                sigcorr_dg_area{iS,kClst,2,i} = sig_corr_dg(in & inct & area(:,1)~=area(:,2));
            end
            kClst = kClst+1;
        end
    end
end

%% plot
clr = {[1 0.6 0.6],[1 0 0];[0.6 0.6 0.6],[0 0 0]};
xticks = {'Nomod vs. Nomod';'Nomod vs. Iact';'Nomod vs. Dact';'Nomod vs. Inh';...
    'Iact vs. Iact';'Iact vs. Dact';'Iact vs. Inh';'Dact vs. Dact';'Dact vs. Inh';'Inh vs. Inh'};
titlelist = {'total';'rs vs. rs';'fs vs. fs';'rs vs. fs'};

for iCT = 1:4
    subplot(1,4,iCT)
    data = cell(2,1);
    hold on;
    for i = 1:2
        data{2/i} = cellfun(@(x) nanmean(x),squeeze(noisecorr_sp_area(:,:,2/i,iCT)));
        plot(1:10,data{2/i},'Color',clr{2/i,1});
    end
    for i = 1:2
        errorbar(1:10,nanmean(data{2/i}),nanstd(data{2/i})./sqrt(sum(~isnan(data{2/i}))),...
            'Color',clr{2/i,2})
    end
    for i = 1:size(noisecorr_sp,2)
        h = double(ttest(data{1}(:,i),data{2}(:,i)));
        if h==1
            text(i,h*0.5,'*','Color','r');
        end
    end
    plot([0 11],[0 0],'k:');
    xlim([0 11]);
    title(titlelist{iCT});
    if iCT==1
        ylabel('Noise correlation');
    end
    set(gca,'XTick',1:10,'XTickLabel',xticks,'XTickLabelRotation',45);
end

%%
msquare = cell(2,2);
x = [1 1 1 1 2 2 2 3 3 4];
y = [1 2 3 4 2 3 4 3 4 4];
for iData = 1:2
    for iArea = 1:2
        switch iData
            case 1
                data = cellfun(@nanmean,sigcorr_dg_area(:,:,iArea,1));
                out = sum(cellfun(@(x) size(x,1)<5,sigcorr_dg_area(:,:,iArea,1)),2)>0;
            case 2
                data = cellfun(@nanmean,noisecorr_sp_area(:,:,iArea,1));
                out = sum(cellfun(@(x) size(x,1)<5,noisecorr_sp_area(:,:,iArea,1)),2)>0;
        end
        m = nanmean(data(~out,:));
        msquare{iData,iArea} = NaN(4,4);
        for i = 1:10
            msquare{iData,iArea}(x(i),y(i)) = m(i);
            msquare{iData,iArea}(y(i),x(i)) = m(i);
        end
    end
end

%%
dataList = {'Gratings-signal';'Spontaneous'};
areaList = {'within-area';'between-areas'};
c = cbrewer('div','RdBu',50);
c([1:5, 46:50],:) = [];
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 6]);
for iData = 1:2
    for iArea = 1:2
        axes('Position',axpt(2,2,iArea,iData,axpt(15,10,2:15,1:9)));
        imagesc(msquare{iData,iArea});
        colormap(c);
        axis xy
        set(gca,'CLim',[-0.15 0.15],'FontSize',5,'LineWidth',0.35,'XTick',1:4,...
            'XTickLabel',{'Nomod';'I act';'D act';'Inh'},'YTick',1:4,'YTickLabel',...
            {'Nomod';'I act';'D act';'Inh'},'XTickLabelRotation',45,'Box','off','TickDir','out');
        if iData==1
            set(gca,'XTickLabel',[]);
            title(areaList{iArea});
        end
        if iArea==1
            ylabel(dataList{iData},'FontSize',5);
        else
            set(gca,'YTickLabel',[]);
        end
    end
end
print(fHandle,'-dtiff','-r600',...
    ['D:\OneDrive - University of California, San Francisco\figures\allen\correlation\correlation_rsrs_',...
    num2str(binsize*1000),'bin_wo_ripple_areas.tif'])




%% subfunctions
function [pairid,noisecorr,time,pupil] = noisecorr_nostim(spikeTime,win,binsize,unit_id,rippletime,rippleon,pupiltime,pupilsize)

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
switch rippleon
    case 1 % w/o ripple events
        spkhist(:,ripplehist>0) = [];
        pupilhist(ripplehist>0) = [];
    case 2 % w/ ripple events
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



function sessionDir = sdir(sessionid)

sessionDir = ['D:\OneDrive\1.allen-andermann\session_',num2str(sessionid),...
    '\session_',num2str(sessionid)];
end