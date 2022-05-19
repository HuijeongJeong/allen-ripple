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
crosscorr = cell((nclst+1+2+1)*(nclst+1+2)/2 + nclst,3,3);
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
    indca1 = ismember(T.unit_id,tag.info.unit_id(strcmp(tag.info.structure,'CA1') &...
        tag.info.ccf(:,3)<7800 & tag.celltype.rs));
    inica1 = ismember(T.unit_id,tag.info.unit_id(strcmp(tag.info.structure,'CA1') &...
        tag.info.ccf(:,3)>8800 & tag.celltype.rs));
    
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
    in = invis|inthal|indca1|inica1;
    spkhist = NaN(length(T.unit_id),nbin);
    spkhist(in,:) =...
        cell2mat(cellfun(@(x) histcounts(x,spontaneous_win(1)-2*win_corr:binsize:...
        spontaneous_win(2)+2*win_corr)*(1/binsize),T.spike_time(in),...
        'UniformOutput',false));
    spkhistz = zscore(spkhist,[],2);
    
    data = NaN(nclst*2+3+1,nbin);
    datasub = NaN(nclst*2,2,nbin);
    datact = NaN(nclst*2+3+1,1);
    kClst = 1;
    for iClst = 1:nclst+3
        for iCT = 1:2
            if iClst>=nclst+1
                switch iClst
                    case nclst+1
                        inclst = inthal;
                    case nclst+2
                        inclst = indca1;
                    case nclst+3
                        inclst = inica1;
                end
                
                if iCT==2
                    continue;
                end
            else
                inclst = invis & celltype(:,iCT) & clstidx==iClst;
                datact(kClst) = iCT;
            end
            data(kClst,:) = conv2(nanmean(spkhistz(inclst,:),1),...
                fspecial('Gaussian',[1 5*resolution],resolution),'same');
            if iClst<nclst+1
                idxtemp = randsample(find(inclst),sum(inclst));
                datasub(iCT+(iClst-1)*2,1,:) = conv(nanmean(spkhistz(idxtemp(1:round(length(idxtemp)/2)),:),1),...
                    fspecial('Gaussian',[1 5*resolution],resolution),'same');
                datasub(iCT+(iClst-1)*2,2,:) = conv(nanmean(spkhistz(idxtemp(round(length(idxtemp)/2)+1:end),:),1),...
                    fspecial('Gaussian',[1 5*resolution],resolution),'same');
            end
            kClst = kClst+1;
        end
    end
    
%     %% lfp power
%     for iRef = 1:2
%         if iscell(filtered_lfp.time{iRef})
%             temp_start = cellfun(@(x) find(x>spontaneous_win(1),1,'first'),...
%                 filtered_lfp.time{iRef},'UniformOutput',false);
%             temp_start = find(cellfun(@(x) ~isempty(x),temp_start),1,'first');
%             temp_end = cellfun(@(x) find(x<spontaneous_win(end),1,'last'),...
%                 filtered_lfp.time{iRef},'UniformOutput',false);
%             temp_end = find(cellfun(@(x) ~isempty(x),temp_end),1,'last');
%             lfptime = cell2mat(filtered_lfp.time{iRef}(temp_start:temp_end)')';
%             flfp = nanmean(cell2mat(filtered_lfp.lfp{iRef}(temp_start:temp_end)),2);
%         else
%             lfptime = filtered_lfp.time{iRef};
%             flfp = nanmean(filtered_lfp.lfp{iRef},2);
%         end
%         
%         data(2*nclst+1+iRef,:) = conv2(interp1(lfptime,abs(flfp)',spktime),...
%             fspecial('Gaussian',[1 5*resolution],resolution),'same');
%     end
%     
    %% pupil
    pupil_size = pupil_data.pupil_height/2.*pupil_data.pupil_width/2*pi;
    pupil_size = interp1(pupil_data.time,pupil_size,spktime);
    pupil_size_norm = pupil_size/nanmean(pupil_size(~inspont));
    data(2*nclst+3+1,:) =...
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
        
        for iPair = 1:3 % rs-rs, fs-fs, rs-fs
            if iPair<3
                [datasub1,datasub2] = deal(data(datact==iPair | isnan(datact),:));
            else
                datasub1 = data(datact==1 | isnan(datact),:);
                datasub2 = data(datact==2 | isnan(datact),:);
            end
            
            k = 1;
            for i = 1:nclst+1+2
                for j = i:nclst+1+2+1
                    if i==j & i>nclst
                        continue;
                    end
                    if iS==1
                        crosscorr{k,iState,iPair} = NaN(nS,2*win_corr/binsize+1);
                        if iState==1 & iPair==1
                            pair = [pair;i,j];
                        end
                    end
                    
                    
                    if iPair==3 & (i>nclst | j>nclst)
                        k = k+1;
                        continue;
                    end
                    corrtemp = NaN(size(analwin,1),2*win_corr/binsize+1);
                    for iw = 1:size(analwin,1)
                        if iPair<3 & (i==j)
                            datatemp = squeeze([datasub(iPair+(i-1)*2,1,analwin(iw,1):analwin(iw,2));...
                                datasub(iPair+(j-1)*2,2,analwin(iw,1):analwin(iw,2))]);
                        else
                            datatemp = [datasub1(i,analwin(iw,1):analwin(iw,2));...
                                datasub2(j,analwin(iw,1):analwin(iw,2))];
                        end
                        idx = 1/binsize+1:size(datatemp,2)-1/binsize;
                        corrtemp(iw,:) = cellfun(@(x) corr(datatemp(1,idx)',datatemp(2,idx+x)'),...
                            num2cell(-win_corr/binsize:win_corr/binsize));
                    end
                    crosscorr{k,iState,iPair}(iS,:) = nanmean(corrtemp);
                    k = k+1;
                end
            end
        end
    end
    
end

figname = {'rsrs';'fsfs';'rsfs'};
ct = {[0.8 0.8 0.8],[0 0 0];[1 0.8 0.8],[1 0 0];[0.8 0.8 1],[0 0 1]};
dataList = {'Iact';'Dact';'Inh';'Thal';'dCA1';'iCA1';'Pupil'};
for iCT = 1:3
    data = squeeze(crosscorr(:,:,iCT));
    fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 10 10]);
    for iP = 1:size(pair,1)
        if isempty(data{iP,iState})
            continue;
        end
        axes('Position',axpt(nclst+2+2,nclst+2+2,pair(iP,1),pair(iP,2)));
        hold on;
        for iState = 1:3
            m = nanmean(data{iP,iState});
            s = nanstd(data{iP,iState})/sqrt(nS);
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
    print(fHandle,'-dtiff','-r600',['cross_corr_around_ripple_',figname{iCT},'.tif']);
end


iCT = 1;
iState = 1;
data = squeeze(crosscorr(:,:,iCT));
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 5]);
for iP = 1:size(pair,1)
    if isempty(data{iP,iState}) | sum(pair(iP,:)>3)>0
        continue;
    end
    axes('Position',axpt(nclst,nclst,pair(iP,1),pair(iP,2),axpt(10,1,2:10,1)));
    hold on;
    m = nanmean(data{iP,iState});
    s = nanstd(data{iP,iState})/sqrt(nS);
    plot(-win_corr:binsize:win_corr,data{iP,iState},'Color',ct{iState,1},'LineWidth',0.35);
%     fill([-win_corr:binsize:win_corr,flip(-win_corr:binsize:win_corr)],...
%         [m+s flip(m-s)],ct{iState,2},'EdgeColor','none');
    plot(-win_corr:binsize:win_corr,m,'Color',ct{iState,2});
    plot([0 0],[-0.8 0.8],'k:');
    ylim([-0.4 0.7])
    if pair(iP,1)==1
        ylabel({['Y: ',dataList{pair(iP,2)}];'corr (X->Y)'});
    else
        set(gca,'YTick',[]);
    end
    if pair(iP,2)==3
        xlabel(['X: ',dataList{pair(iP,1)}]);
    else
        set(gca,'XTick',[]);
    end
    set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35);
    alpha(0.2);
end

cd('D:\OneDrive - University of California, San Francisco\figures\allen\correlation');
print(fHandle,'-dtiff','-r600','cross_corr_rsrs_woripple.tif');


iCT = 1;
iState = 1;
data = squeeze(crosscorr(:,:,iCT));
ii = [1,2,3];
jj = [4,5,6];
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 5]);
for i = 1:3
    for j = 1:3
    axes('Position',axpt(3,3,i,j,axpt(10,1,2:10,1)));
    hold on;
    inpair = pair(:,1)==ii(i) & pair(:,2)==jj(j);
    m = nanmean(data{inpair,iState});
    s = nanstd(data{inpair,iState})/sqrt(nS);
    plot(-win_corr:binsize:win_corr,data{inpair,iState},'Color',ct{iState,1},'LineWidth',0.35);
%     fill([-win_corr:binsize:win_corr,flip(-win_corr:binsize:win_corr)],...
%         [m+s flip(m-s)],ct{iState,2},'EdgeColor','none');
    plot(-win_corr:binsize:win_corr,m,'Color',ct{iState,2});
    plot([0 0],[-0.8 0.8],'k:');
    ylim([-0.4 0.7])
    if i==1
        ylabel({['Y: ',dataList{jj(j)}];'corr (X->Y)'});
    else
        set(gca,'YTick',[]);
    end
    if j==3
        xlabel(['X: ',dataList{ii(i)}]);
    else
        set(gca,'XTick',[]);
    end
    set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35);
    end
end

cd('D:\OneDrive - University of California, San Francisco\figures\allen\correlation');
print(fHandle,'-dtiff','-r600','cross_corr_rs_other_woripple.tif');

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