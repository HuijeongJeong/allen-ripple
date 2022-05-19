clc; clearvars; close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

refList = {'medial';'lateral';'global'};

nclst = 3;
resolution = 2;
win_corr = 1;
binsize = 0.01;
iCT = 2;

[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
session_id = tag.info.session_id(idx);

unitid = unit_id.vis;
sessionList = unique(session_id);
nS = length(sessionList);

pair = [];
crosscorr = cell((nclst+1+2+1)*(nclst+1+2)/2,3,3);
for iS = 1:nS
    iS
    %% load data
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_ripples.mat'],'CA1_ripple_classified',...
        'filtered_lfp','spontaneous_anal_win','spontaneous_win');
    load([sdir(sessionList(iS)),'_Events.mat'],'running_speed','pupil_data');
    
    %% calculate index
    invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis));
    inthal = ismember(T.unit_id,tag.info.unit_id(ismember(tag.info.structure,{'LGd';'LGv';'LP'})));
    
    [in,idx] = ismember(T.unit_id,unitid);
    clstidx = NaN(length(T.unit_id),1);
    clstidx(in) = cluster_idx.vis{nclst-1}(idx(in));
    
    [~,idx] = ismember(T.unit_id,tag.info.unit_id);
    celltype = [tag.celltype.rs(idx), tag.celltype.fs(idx)];
    
    %% analysis window
    spktime = spontaneous_win(1)-2*win_corr+binsize/2:binsize:spontaneous_win(2)+2*win_corr-binsize/2;
    nbin = length(spktime);
    
    inripple = false(length(spktime),3);
    for iRef = 1:3
        for iRp = 1:size(CA1_ripple_classified.(refList{iRef}),1)
            [~,idx] = min(abs(spktime-CA1_ripple_classified.(refList{iRef})(iRp,1)));
            inripple(idx-1.5/binsize:idx+1.5/binsize,iRef) = true;
        end
    end
    
    inspont = false(length(spktime),1);
    for i = 1:size(spontaneous_anal_win,1)
        inspont(spktime>=spontaneous_anal_win(i,1) & spktime<=spontaneous_anal_win(i,2)) = true;
    end
    
    %% average psth for each cluster
    spkhist = NaN(length(T.unit_id),nbin);
    spkhist(invis|inthal,:) =...
        cell2mat(cellfun(@(x) histcounts(x,spontaneous_win(1)-2*win_corr:binsize:...
        spontaneous_win(2)+2*win_corr)*(1/binsize),T.spike_time(invis|inthal),...
        'UniformOutput',false));
    spkhistz = zscore(spkhist,[],2);
    
    data = NaN(nclst*2+1+2+1,nbin);
    datact = NaN(nclst*2+1+2+1,1);
    kClst = 1;
    for iClst = 1:nclst+1
        for iCT = 1:2
            if iClst==1
                inclst = inthal;
                if iCT==2
                    continue;
                end
            else
                inclst = invis & celltype(:,iCT) & clstidx==iClst-1;
                datact(kClst) = iCT;
            end
            data(kClst,:) = conv2(nanmean(spkhistz(inclst,:),1),...
                fspecial('Gaussian',[1 5*resolution],resolution),'same');
            kClst = kClst+1;
        end
    end
    
    %% lfp power
    for iRef = 1:2
        if iscell(filtered_lfp.time{iRef})
            temp_start = cellfun(@(x) find(x>spontaneous_win(1),1,'first'),...
                filtered_lfp.time{iRef},'UniformOutput',false);
            temp_start = find(cellfun(@(x) ~isempty(x),temp_start),1,'first');
            temp_end = cellfun(@(x) find(x<spontaneous_win(end),1,'last'),...
                filtered_lfp.time{iRef},'UniformOutput',false);
            temp_end = find(cellfun(@(x) ~isempty(x),temp_end),1,'last');
            lfptime = cell2mat(filtered_lfp.time{iRef}(temp_start:temp_end)')';
            flfp = nanmean(cell2mat(filtered_lfp.lfp{iRef}(temp_start:temp_end)),2);
        else
            lfptime = filtered_lfp.time{iRef};
            flfp = nanmean(filtered_lfp.lfp{iRef},2);
        end
        
        data(2*nclst+1+iRef,:) = conv2(interp1(lfptime,abs(flfp)',spktime),...
            fspecial('Gaussian',[1 5*resolution],resolution),'same');
    end
    
    %% pupil 
    pupil_size = pupil_data.pupil_height/2.*pupil_data.pupil_width/2*pi;
    pupil_size = interp1(pupil_data.time,pupil_size,spktime);
    pupil_size_norm = pupil_size/nanmean(pupil_size(~inspont));
    data(2*nclst+1+2+1,:) =...
        conv(pupil_size_norm,fspecial('Gaussian',[1 5*resolution],resolution),'same');
    
    %% calculating correlation
    for iState = 1:3 %immobile w/o ripple, immobile w/ dorsal ripple, w/ intermediate ripple
        if iState==1
            inanal = inspont & sum(inripple,2)==0;
        else
            inanal = inspont & inripple(:,iState-1);
        end
        analwin = [find(diff([0;inanal])==1),find(diff([inanal;0])==-1)];
        analwin(diff(analwin,[],2)<3/binsize,:) = []; % use only window that is longer than 3s
        
        for iPair = 1:2 % rs-rs, fs-fs
            datasub = data(datact==iPair | isnan(datact),:);
            
            k = 1;
            for i = 1:nclst+1+2
                for j = i+1:nclst+1+2+1
                    corrtemp = NaN(size(analwin,1),2*win_corr/binsize+1);
                    for iw = 1:size(analwin,1)
                        datatemp = [datasub(i,analwin(iw,1):analwin(iw,2));...
                            datasub(j,analwin(iw,1):analwin(iw,2))];
                        idx = 1/binsize+1:size(datatemp,2)-1/binsize;                        
                        corrtemp(iw,:) = cellfun(@(x) corr(datatemp(1,idx)',datatemp(2,idx+x)'),...
                            num2cell(-win_corr/binsize:win_corr/binsize));
                    end
                    if iS==1
                        crosscorr{k,iState,iPair} = NaN(nS,2*win_corr/binsize+1);
                        if iState==1 & iPair==1
                            pair = [pair;i,j];
                        end
                    end
                    crosscorr{k,iState,iPair}(iS,:) = nanmean(corrtemp);
                    k = k+1;
                end
            end
        end
    end
    
    
     %% calculating correlation
    for iState = 1:3 %immobile w/o ripple, immobile w/ dorsal ripple, w/ intermediate ripple
        if iState==1
            inanal = inspont & sum(inripple,2)==0;
        else
            inanal = inspont & inripple(:,iState-1);
        end
        analwin = [find(diff([0;inanal])==1),find(diff([inanal;0])==-1)];
        analwin(diff(analwin,[],2)<3/binsize,:) = []; % use only window that is longer than 3s
        
        for iPair = 1:2 % rs-rs, fs-fs
            datasub = data(datact==iPair | isnan(datact),:);
            
            k = 1;
            for i = 1:nclst+1+2
                for j = i+1:nclst+1+2+1
                    corrtemp = NaN(size(analwin,1),2*win_corr/binsize+1);
                    for iw = 1:size(analwin,1)
                        datatemp = [datasub(i,analwin(iw,1):analwin(iw,2));...
                            datasub(j,analwin(iw,1):analwin(iw,2))];
                        idx = 1/binsize+1:size(datatemp,2)-1/binsize;                        
                        corrtemp(iw,:) = cellfun(@(x) corr(datatemp(1,idx)',datatemp(2,idx+x)'),...
                            num2cell(-win_corr/binsize:win_corr/binsize));
                    end
                    if iS==1
                        crosscorr{k,iState,iPair} = NaN(nS,2*win_corr/binsize+1);
                        if iState==1 & iPair==1
                            pair = [pair;i,j];
                        end
                    end
                    crosscorr{k,iState,iPair}(iS,:) = nanmean(corrtemp);
                    k = k+1;
                end
            end
        end
    end
    
    
end



ct = {[0.8 0.8 0.8],[0 0 0];[1 0.8 0.8],[1 0 0];[0.8 0.8 1],[0 0 1]};
dataList = {'Thal';'Iact';'Dact';'Inh';'D Rp power';'I Rp power';'Pupil'};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 10 10]);
for iP = 1:size(pair,1)
    axes('Position',axpt(nclst+1+2,nclst+1+2,pair(iP,1),pair(iP,2)-1));
    hold on;
    
    for iState = 1:3
        m = nanmean(crosscorr{iP,iState});
        s = nanstd(crosscorr{iP,iState})/sqrt(nS);
        fill([-win_corr:binsize:win_corr,flip(-win_corr:binsize:win_corr)],...
            [m+s flip(m-s)],ct{iState,2},'EdgeColor','none');
        plot(-win_corr:binsize:win_corr,m,'Color',ct{iState,2});
    end
    plot([0 0],[-0.8 0.8],'k:');
    ylim([-0.3 0.7])
    if pair(iP,1)==1
        ylabel({['Y: ',dataList{pair(iP,2)}];'corr (X->Y)'});
    else
        set(gca,'YTick',[]);
    end
    if pair(iP,2)==nclst+1+2+1
        xlabel(['X: ',dataList{pair(iP,1)}]);
    else
        set(gca,'XTick',[]);
    end
    set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35);
    alpha(0.2);
end

cd('D:\OneDrive - University of California, San Francisco\figures\allen\correlation');
print(fHandle,'-dtiff','-r600','cross_corr_around_ripple_rsrs.tif');

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