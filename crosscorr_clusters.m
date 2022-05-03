clc; clearvars; close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

refList = {'medial';'lateral'};
typeList = {'rs';'fs'};

nclst = 3;
resolution = 10;
win = [-3 3];
win_corr = [-1 1];
binsize = 0.01;
time = win(1)+binsize/2:binsize:win(2)-binsize/2;

[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
session_id = tag.info.session_id(idx);

unitid = unit_id.vis;
sessionList = unique(session_id);
nS = length(sessionList);

pair = [];
crosscorr = cell(10,2);
for iS = 19:nS
    iS
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_ripple_modulated.mat'],'fr_ripple');
    load([sdir(sessionList(iS)),'_ripples.mat'],'CA1_ripple_classified','filtered_lfp');
    load([sdir(sessionList(iS)),'_Events.mat'],'running_speed','pupil_data');
    
    invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis));
    inthal = ismember(T.unit_id,tag.info.unit_id(ismember(tag.info.structure,{'LGd';'LGv';'LP'})));
    
    [in,idx] = ismember(T.unit_id,unitid);
    clstidx = NaN(length(T.unit_id),1);
    clstidx(in) = cluster_idx.vis{nclst-1}(idx(in));
    
    [~,idx] = ismember(T.unit_id,tag.info.unit_id);
    celltype = [tag.celltype.rs(idx), tag.celltype.fs(idx)];
    
    if ismember('hist',fieldnames(fr_ripple))
        have_fr = 1;
    else
        have_fr = 0;
    end
    
    
    for iRef = 1:2
        data = cell(5,1);
        rippletime = CA1_ripple_classified.(refList{iRef})(:,1);
        if have_fr==1
            intime = fr_ripple.time>=win(1) & fr_ripple.time<=win(2);
            spkhist = fr_ripple.hist{iRef};
            spkhist = cellfun(@(x) x(:,intime),spkhist,'UniformOutput',false);
        else
            [spkhist,spikeTimes] = deal(cell(length(T.unit_id),1));
            spikeTimes(invis|inthal) = cellfun(@(x) spikeWin(x,rippletime,win),...
                T.spike_time(invis|inthal),'UniformOutput',false);
            [~,spkhist(invis|inthal)] = cellfun(@(x) spikeBin(x,win,binsize,binsize),...
                spikeTimes(invis|inthal),'UniformOutput',false);
        end
        frconv = cellfun(@(x) conv2(x,fspecial('Gaussian',[1 5*resolution],resolution),'same'),...
            spkhist,'UniformOutput',false);
        frconvz = cellfun(@(x) (x-mean(x(:)))/std(x(:)),frconv,'UniformOutput',false);
        
        data{1} = squeeze(mean(cell2mat(permute(frconvz(inthal),[3,2,1])),3));
        for iclst = 1:3
            in = invis & celltype(:,1) & clstidx==iclst;
            data{iclst+1} = squeeze(mean(cell2mat(permute(frconvz(in),[3,2,1])),3));
        end
        if iscell(filtered_lfp.time{iRef})
            temp_start = cellfun(@(x) find(x>rippletime(1),1,'first'),...
                filtered_lfp.time{iRef},'UniformOutput',false);
            temp_start = find(cellfun(@(x) ~isempty(x),temp_start),1,'first');
            temp_end = cellfun(@(x) find(x<rippletime(end),1,'last'),...
                filtered_lfp.time{iRef},'UniformOutput',false);
            temp_end = find(cellfun(@(x) ~isempty(x),temp_end),1,'last');
            lfptime = cell2mat(filtered_lfp.time{iRef}(temp_start:temp_end)')';
            flfp = nanmean(cell2mat(filtered_lfp.lfp{iRef}(temp_start:temp_end)),2);
        else
            lfptime = filtered_lfp.time{iRef};
            flfp = nanmean(filtered_lfp.lfp{iRef},2);
        end
        [time_lfp,lfp] = alignbeh2event(lfptime,flfp,rippletime,win);
        data{5} = interp1(time_lfp,abs(lfp)',time)';
        
        pupildata = conv(pupil_data.pupil_height/2.*pupil_data.pupil_width/2*pi,...
            fspecial('Gaussian',[1 5*3],3),'same');
        [time_pupil,pupil] = alignbeh2event(pupil_data.time,pupildata,rippletime,win);
        data{6} = interp1(time_pupil,pupil',time)';
        
        k = 1;
        for i = 1:5
            for j = i+1:6
                idx = find(time>=-1.5 & time<=1.5);
                if iS==1
                    crosscorr{k,iRef} = NaN(nS,diff(win_corr)/binsize+1);
                end
                if ~isempty(data{i}) & ~isempty(data{j})
                    inpair = sum(isnan(data{i}),2)==0 & sum(isnan(data{j}),2)==0;
                    crosscorr{k,iRef}(iS,:) =...
                        nanmean(cell2mat(cellfun(@(y) cellfun(@(x) corr(data{i}(y,idx)',data{j}(y,idx+x)'),...
                        num2cell(win_corr(1)/binsize:win_corr(2)/binsize)),...
                        num2cell(find(inpair)),'UniformOutput',false)));
                end
                if iS==1 & iRef==1
                    pair = [pair;i,j];
                end
                k = k+1;
            end
        end
    end
end


ct = {[1 0.8 0.8],[1 0 0];[0.8 0.8 1],[0 0 1]};
dataList = {'Thal';'Iact';'Dact';'Inh';'Rp power'};
figure;
for iP = 1:size(pair,1)
    axes('Position',axpt(4,4,pair(iP,1),pair(iP,2)-1));
    hold on;
    for iRef = 1:2
        plot(win_corr(1):binsize:win_corr(2),crosscorr{iP,iRef},'Color',ct{iRef,1},'LineWidth',0.35);
    end
    for iRef = 1:2
        plot(win_corr(1):binsize:win_corr(2),nanmean(crosscorr{iP,iRef}),'Color',ct{iRef,2},'LineWidth',1);
    end
    plot([0 0],[-0.8 0.8],'k:');
    ylim([-0.7 0.7])
    if pair(iP,1)==1
       ylabel({['Y: ',dataList{pair(iP,2)}];'cross-corr (X->Y)'});
    else
        set(gca,'YTick',[]);
    end
    if pair(iP,2)==5
       xlabel(['X: ',dataList{pair(iP,1)}]);
    else
        set(gca,'XTick',[]);        
    end
end

function [time,behavior] = alignbeh2event(timeline,behdata,eventtime,win)
[~,minidx] = cellfun(@(x) min(abs(timeline-x)),num2cell(eventtime));
binsize = nanmean(diff(timeline));
winidx = [floor(win(1)/binsize) ceil(win(2)/binsize)];
time = [winidx(1):winidx(2)]*binsize;
behavior = cell2mat(cellfun(@(x) behdata([winidx(1):winidx(2)]+x)',num2cell(minidx),'UniformOutput',false));
end


function sessionDir = sdir(sessionid)

sessionDir = ['D:\OneDrive\1.allen-andermann\session_',num2str(sessionid),...
    '\session_',num2str(sessionid)];
end