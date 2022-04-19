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
noisecorr_sp = cell(nS,k*(k+1)/2+k+1,4); 

% noise correlation during spontaneous period within each area or across areas
% (session, cluster pair, within/across areas, celltype)
noisecorr_sp_area = cell(nS,k*(k+1)/2+k+1,2,4); 

binsize = 2; % unit: second

%% 
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_Events.mat'],'running_speed','pupil_data');
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_anal_win','spontaneous_CA1_ripple');
    
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

function sessionDir = sdir(sessionid)

sessionDir = ['D:\OneDrive\1.allen-andermann\session_',num2str(sessionid),...
    '\session_',num2str(sessionid)];
end