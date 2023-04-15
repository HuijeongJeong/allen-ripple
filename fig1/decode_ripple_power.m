% clearvars; clc; close all;

rng(3);

load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');

%%
nclst =3;
binsize = 0.2;
niter = 100;
t = 11;

unitid = tag.info.unit_id(tag.area.vis);
celltype = [tag.celltype.rs(tag.area.vis),tag.celltype.fs(tag.area.vis)];
[in,idx] = ismember(unitid,unit_id.vis);
clusteridx = zeros(length(unitid),1);
clusteridx(in) = cluster_idx.vis{nclst-1}(idx(in));

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

rpType = {'medial','lateral','global'};

% [accuracy,accuracy_sf] = deal(cell(3,2));
% accuracy_beh = nan(nS,t);
% nCell = cell(2,1);
% nRipple = nan(nS,1);
%%
for iS = 19:nS
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_Events'],'running_speed','pupil_data');
    load([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_win','spontaneous_anal_win',...
        'CA1_ripple_classified','filtered_lfp');

    %%
    [in,idx] = ismember(T.unit_id,unitid);
    spktime = T.spike_time(in);
    clstidx = clusteridx(idx(in));
    ctype = celltype(idx(in),:);

    spkbin = cell2mat(cellfun(@(x) histcounts(x,spontaneous_win(1):binsize:spontaneous_win(2)),...
        spktime,'UniformOutput',false));

    [~,~,pupilidx] = histcounts(pupil_data.time,spontaneous_win(1):binsize:spontaneous_win(2));
    pupilsize = pupil_data.pupil_width/2.*pupil_data.pupil_height/2*pi;
    pupilbin = cellfun(@(x) mean(pupilsize(pupilidx==x)),num2cell(1:max(pupilidx)));

    [~,~,runningidx] = histcounts(running_speed.time,spontaneous_win(1):binsize:spontaneous_win(2));
    runningbin = cellfun(@(x) mean(running_speed.velocity(runningidx==x)),num2cell(1:max(runningidx)));

    %%
    nssbin = NaN(floor(diff(spontaneous_win)/binsize),2);
    for iRp = 1:2
        if ~iscell(filtered_lfp.lfp{iRp})
            nss = mean(NormalizedSquaredSignal_HJ(...
                [filtered_lfp.time{iRp}',filtered_lfp.lfp{iRp}]),2);
            [~,~,nssidx] = histcounts(filtered_lfp.time{iRp},spontaneous_win(1):binsize:spontaneous_win(2));
        else
            nss = cell(length(filtered_lfp.lfp{iRp}),1);
            for iW = 1:length(filtered_lfp.lfp{iRp})
                nss{iW} = mean(NormalizedSquaredSignal_HJ(...
                    [filtered_lfp.time{iRp}{iW}',filtered_lfp.lfp{iRp}{iW}]),2);
            end
            nss = cell2mat(nss);
            [~,~,nssidx] = histcounts(cell2mat(filtered_lfp.time{iRp}'),spontaneous_win(1):binsize:spontaneous_win(2));
        end
        nssbin(:,iRp) = cellfun(@(x) mean(nss(nssidx==x)),num2cell(1:max(nssidx)));
    end

    %%
    ripplebin = zeros(floor(diff(spontaneous_win)/binsize),3);
    for iRp = 1:3
        [~,~,rpidx] = histcounts(CA1_ripple_classified.(rpType{iRp})(:,1),...
            spontaneous_win(1):binsize:spontaneous_win(2));
        ripplebin(unique(rpidx),iRp) = 1;
    end

    %%
    rippleidx = find(sum(ripplebin,2)>0 & movmean(double(sum(ripplebin,2)>0),t,'Endpoints','fill')==1/t);
    nonrippleidx = find(sum(ripplebin,2)==0 & movmean(double(sum(ripplebin,2)>0),t,'Endpoints','fill')==0);
    nripple = length(rippleidx);
    

    [score,score_sf] = deal(cell(3,2));
    score_beh = nan(niter,t);

    for iiter = 1:niter
        fprintf('session #%d, iteration #%d\n',iS,iiter);
        iteridx = [randsample(rippleidx,nripple),randsample(nonrippleidx,nripple)];
        trainidx = iteridx(1:round(nripple/2),:);
        testidx = iteridx(round(nripple/2)+1:end,:);
        trainrsp = repmat([1,2],round(nripple/2),1);
        testrsp = repmat([1,2],nripple-round(nripple/2),1);

        for iCT = 1:3
            switch iCT
                case 1
                    incelltype = true(size(spkbin,1),1);
                case {2,3}
                    incelltype = ctype(:,iCT-1);
            end
            for i = 1:2
                if isempty(score{iCT,i})
                    score{iCT,i} = NaN(niter,t);
                    score_sf{iCT,i} = NaN(niter,t);
                end
                if i==1
                    spktemp = spkbin(incelltype & clstidx==0,:);
                else
                    spktemp = spkbin(incelltype & clstidx>0,:);
                end

                for it = 1:t
                    spktrain = spktemp(:,trainidx(:)+(it-round(t/2)));
                    mdl = fitcsvm(spktrain',trainrsp(:));
                    mdlsf = fitcsvm(spktrain',randsample(trainrsp(:),length(trainrsp(:))));
                    score{iCT,i}(iiter,it) = mean(predict(mdl,spktemp(:,testidx(:)+(it-round(t/2)))')==testrsp(:));
                    score_sf{iCT,i}(iiter,it) = mean(predict(mdlsf,spktemp(:,testidx(:)+(it-round(t/2)))')==...
                        randsample(testrsp(:),length(testrsp(:))));

                    if iCT==3 & i==2
                        behtrain = [runningbin(trainidx(:)+(it-round(t/2)))',pupilbin(trainidx(:)+(it-round(t/2)))'];
                        behtest = [runningbin(testidx(:)+(it-round(t/2)))',pupilbin(testidx(:)+(it-round(t/2)))'];
                        mdl = fitcsvm(behtrain,trainrsp(:));
                        score_beh(iiter,it) = mean(predict(mdl,behtest)==testrsp(:));
                    end
                end
            end

        end
    end
    accuracy = cellfun(@(x,y) [x;nanmean(y)],accuracy,score,'UniformOutput',false);
    accuracy_sf = cellfun(@(x,y) [x;nanmean(y)],accuracy_sf,score_sf,'UniformOutput',false);
    accuracy_beh(iS,:) = mean(score_beh);
    nRipple(iS) = nripple;
    for iCT = 1:2
        if isempty(nCell{iCT})
            nCell{iCT} = nan(nS,4);
        end
        nCell{iCT}(iS,:) = sum([clstidx==0 clstidx==1 clstidx==2 clstidx==3] & ctype(:,iCT));
    end
end

%%
ct = {'k','r'};
celltypeList = {'RS+FS';'RS';'FS'};

figure
for iCT = 1:3
    if iCT==1
out = sum(nCell{1}+nCell{2}==0,2)>0;
    else
out = sum(nCell{iCT-1}==0,2)>0;
    end
subplot(1,4,iCT)
hold on;
for i = 1:2
    fill([-1:0.2:1,flip(-1:0.2:1)],[mean(accuracy{iCT,i}(~out,:))+std(accuracy{iCT,i}(~out,:))/sqrt(sum(~out)),...
        flip(mean(accuracy{iCT,i}(~out,:))-std(accuracy{iCT,i}(~out,:))/sqrt(sum(~out)))],...
        ct{i},'EdgeColor','none')
    fill([-1:0.2:1,flip(-1:0.2:1)],[mean(accuracy_sf{iCT,i}(~out,:))+std(accuracy_sf{iCT,i}(~out,:))/sqrt(sum(~out)),...
        flip(mean(accuracy_sf{iCT,i}(~out,:))-std(accuracy_sf{iCT,i}(~out,:))/sqrt(sum(~out)))],...
        ct{i},'EdgeColor','none')
plot([-1:0.2:1],mean(accuracy{iCT,i}(~out,:)),'Color',ct{i})
plot([-1:0.2:1],mean(accuracy_sf{iCT,i}(~out,:)),'Color',ct{i},'LineStyle',':')
end
[~,p] = ttest(accuracy{iCT,1}(~out,:),accuracy{iCT,2}(~out,:));
h = double(p<0.05);
h(h==0) = NaN;
plot(-1:0.2:1,h*0.78,'k','LineWidth',2)
plot([0 0],[0.45 0.8],'k--')
ylim([0.45 0.8]);
alpha(0.2)
title(celltypeList{iCT})
xlabel('lag (s)')
end

subplot(1,4,4)
hold on;
fill([-1:0.2:1 flip(-1:0.2:1)],[mean(accuracy_beh(~out,:))+std(accuracy_beh(~out,:))/sqrt(sum(~out)),...
    flip(mean(accuracy_beh(~out,:))-std(accuracy_beh(~out,:)/sqrt(sum(~out))))],'k','EdgeColor','none')
plot(-1:0.2:1,mean(accuracy_beh(~out,:)),'k')
plot([0 0],[0.45 0.8],'k--')
ylim([0.45 0.8])
alpha(0.2);
title('Running+pupil')
xlabel('lag (s)')







