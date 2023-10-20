clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');

k = 3;
sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);
celltype = [tag.celltype.rs, tag.celltype.fs];

speedlimit = 2;
binsize = 2;

threshold = [2 5];   %threshold for ripple beginning/end and peak
duration = [30 20 250];    %min inter-ripple inter, min ripple length, max ripple length
Fs = 1250;

[rppower_mean,rppower_peak,rppower_mean_dg,rppower_peak_dg] = deal(cell(nS,4));
[runningave,pupilave,runningave_dg,pupilave_dg] = deal(NaN(nS,4));
%%
typeList = {'rs';'fs'};
iType = 1;
for iS = 1:length(sessionList)
    iS
    clear pupil_data
    load([sdir(sessionList(iS)),'_Events.mat'],'drifting_gratings','running_speed','pupil_data');
    load([sdir(sessionList(iS)),'_ripples.mat'],'filtered_lfp','spontaneous_anal_win','spontaneous_win');
    load([sdir(sessionList(iS)),'_stimulus_evoked.mat'],'beh_gratings');
    
    % detect ripples
    win = [drifting_gratings.window(1,1) drifting_gratings.window(end,end)];
    start = running_speed.time(diff([0;running_speed.immobile])==1);
    stop = running_speed.time(diff(running_speed.immobile)==-1);
    if length(stop)<length(start)
        immobile = [start,[stop;win(2)]];
    else
        immobile = [start,stop];
    end
    
    %%
    [ripples,nss,nsstime] = deal(cell(2,1));
    for iProbe = 1:2
        if iscell(filtered_lfp.time{iProbe})
            ripple_tmp = cell2mat(cellfun(@(x,y) FindRipples_HJ([x',y],...
                'thresholds',threshold,'durations',duration,'frequency',Fs),...
                filtered_lfp.time{iProbe},filtered_lfp.lfp{iProbe},'UniformOutput',false));
            nss_tmp = cell(length(filtered_lfp.time{iProbe}),1);
            for iw = 1:length(filtered_lfp.time{iProbe})
               nss_tmp{iw} =  mean(NormalizedSquaredSignal_HJ([filtered_lfp.time{iProbe}{iw}',...
                filtered_lfp.lfp{iProbe}{iw}]),2);
            end 
            nss{iProbe} = cell2mat(nss_tmp);
            nsstime{iProbe} = cell2mat(filtered_lfp.time{iProbe}');
        else
            ripple_tmp = FindRipples_HJ([filtered_lfp.time{iProbe}',filtered_lfp.lfp{iProbe}],...
                'thresholds',threshold,'durations',duration,'frequency',Fs);
            nss{iProbe} = mean(NormalizedSquaredSignal_HJ([filtered_lfp.time{iProbe}',...
                filtered_lfp.lfp{iProbe}]),2);
            nsstime{iProbe} = filtered_lfp.time{iProbe};
        end
        inImmobile = logical(sum(cell2mat(cellfun(@(x) ripple_tmp(:,1)>x(1) & ripple_tmp(:,3)<x(2),...
            mat2cell(immobile,ones(size(immobile,1),1),2),'UniformOutput',false)'),2));
        ripples{iProbe} = ripple_tmp(inImmobile,:);
    end
    
    %%
    dglist = sort(unique(drifting_gratings.stimulus_condition_id));
    ndg = length(dglist);
    dgIndex = false(length(drifting_gratings.stimulus_condition_id),ndg);
    for idg = 1:ndg
        dgIndex(:,idg) = drifting_gratings.stimulus_condition_id==dglist(idg);
    end

    pupil_size = pupil_data.pupil_height/2.*pupil_data.pupil_width/2*pi;

    ripples = cell2mat(ripples);
    dg_ripple = cellfun(@(x) sum(ripples(:,1)>=x(1) & ripples(:,1)<=x(2)),...
        mat2cell(drifting_gratings.window,ones(size(drifting_gratings.window,1),1),2));
    
    [rppower_mean_dg(iS,:),rppower_peak_dg(iS,:),runningave_dg(iS,:),pupilave_dg(iS,:)] =...
        ripplepower_stim(nss,nsstime,drifting_gratings.window,...
        cell2mat(beh_gratings.velocity_ave'),speedlimit,dg_ripple,...
            cell2mat(beh_gratings.pupil_area_ave'));
        
    [rppower_mean(iS,:),rppower_peak(iS,:),runningave(iS,:),pupilave(iS,:)] =...
        ripplepower_nostim(nss,nsstime,win,binsize,ripples,pupil_data.time,...
        pupil_size,running_speed.time,running_speed.velocity_conv,speedlimit);

end
%%
% [in,idx] = ismember(unit_id.vis,tag.info.unit_id(tag.celltype.rs));
% sessionid = tag.info.session_id(tag.celltype.rs);
% sessionid = sessionid(idx(in));
% cidx = cluster_idx.vis{2}(in);
% n = nan(nS,3);
% for iS = 1:nS
%     for iClst= 1:3
%    n(iS,iClst) = sum(cidx(sessionid==sessionList(iS))==iClst);
%     end
% end
% out = sum(n==0,2)>0;

out = false(length(sessionList),1);
out([7,19]) = true;
%%
dataList = {'spont','dg'};
for iData = 1
    if iData==1
        rave = runningave(~out,:);
        pave = pupilave(~out,:)./repmat(pupilave(~out,4),1,4);
        rp = rppower_peak(~out,:);
    else
        in = sum(cellfun(@isempty,rppower_peak_dg),2)==0;
        rave = runningave_dg(~out&in,:);
        pave = pupilave_dg(~out&in,:)./repmat(pupilave_dg(~out,4),1,4);
        rp = rppower_peak_dg(~out&in,:);
    end
        
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6.5 20/6]);
h(1) = axes('Position',axpt(3,10,1,1:8,[],[0.2 0.05]));
hold on;
plot(1:4,rave,'Color',[0.6 0.6 0.6]);
errorbar(1:4,mean(rave),std(rave)/sqrt(sum(~out)),'k','CapSize',3);
set(gca,'YTick',0:10:30,'YLim',[0 30]);
ylabel('Running speed (cm/s)');
% [tbl,rm] = simple_mixed_anova(rave)
pr = [];
statr = [];
for x = 1:3
    for y = x+1:4
        [~,p,~,stat] = ttest(rave(:,x),rave(:,y));
        pr = [pr;p*6];
        statr = [statr;stat.tstat];
    end
end

h(2) = axes('Position',axpt(3,10,2,1:8,[],[0.2 0.05]));
hold on;
plot(1:4,pave,'Color',[0.6 0.6 0.6]);
errorbar(1:4,mean(pave),std(pave)/sqrt(sum(~out)),'k','CapSize',3);
set(gca,'YTick',0:0.5:1,'YLim',[0 1.2]);
ylabel({'Norm. pupil size';'(1 = during mobile)'});
pp = [];
statp = [];
for x = 1:3
    for y = x+1:4
        [~,p,~,stat] = ttest(pave(:,x),pave(:,y));
        pp = [pp;p*6];
        statp = [statp;stat.tstat];
    end
end

ct = cbrewer('qual','Paired',6);
ct(3:4,:) = [];
ct = flip(ct);
h(3) = axes('Position',axpt(3,10,3,1:8,[],[0.2 0.05]));
hold on;
for iRp = 1:2
    m = cellfun(@(x) nanmean(x(:,iRp)),rp);
   plot(1:4,m,'Color',ct(2*iRp,:));
   errorbar(1:4,mean(m),std(m)/sqrt(size(m,1)),'Color',ct(2*iRp-1,:),'CapSize',3);
   prp = [];
statrp = [];
for x = 1:3
    for y = x+1:4
        [~,p,~,stat] = ttest(m(:,x),m(:,y));
        prp = [prp;p*6];
        statrp = [statrp;stat.tstat];
    end
end

end
set(gca,'YTick',0:5:20,'YLim',[0 22]);
ylabel('Norm. ripple band power');


set(h,'XLim',[0.5 4.5],'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,...
    'XTick',1:4,'XTickLabel',{'Ripple';'Low pupil';'High pupil';'Mobile'},'XTickLabelRotation',45)

cd('D:\OneDrive - University of California, San Francisco\figures\2.allen\fig2');
print(fHandle,'-depsc','-painters',...
    ['ripplepower_',num2str(binsize*1000),'bin_wo_ripple_',dataList{iData},'.ai'])
end

%%
        
function [rppower_mean,rppower_peak,running_ave,pupil_ave] =...
    ripplepower_stim(nss,nsstime,win,speed,speedlimit,trialripple,pupilsize)

[nssmean,nsspeak] = deal(NaN(size(win,1),2));
for iRp = 1:2
    in = cellfun(@(x) nsstime{iRp}>=x(1) & nsstime{iRp}<=x(2),...
        mat2cell(win,ones(size(win,1),1),2),'UniformOutput',false);
    nssmean(:,iRp) = cellfun(@(x) mean(nss{iRp}(x)),in);
    peaktemp = cellfun(@(x) max(nss{iRp}(x)),in,'UniformOutput',false);
    peaktemp(cellfun(@isempty,peaktemp)) = {NaN};
    nsspeak(:,iRp) = cell2mat(peaktemp);
end

[rppower_mean,rppower_peak] = deal(cell(4,1));
[running_ave,pupil_ave] = deal(NaN(4,1));

pupilmobile = nanmean(pupilsize(abs(speed)>speedlimit));
for iState = 1:4
    switch iState
        case 1 % no ripple
            instate = trialripple>0;
        case {2,3}% immobile & no ripple & low (or high) pupil (bottom (or top) 50% pupil within immobile&no ripple condition)
            
            inconditiontemp = abs(speed)<speedlimit & trialripple==0;
            if iState==2
                instate = inconditiontemp & pupilsize<pupilmobile*0.5;
            elseif iState==3
                instate = inconditiontemp & pupilsize>=pupilmobile*0.5;
            end
        case 4 % mobile
            instate = abs(speed)>speedlimit;
    end
    if sum(instate)==0
        continue
    end
    rppower_mean{iState} = nssmean(instate,:);
    rppower_peak{iState} = nsspeak(instate,:);
    running_ave(iState) = mean(speed(instate));
    pupil_ave(iState) = nanmean(pupilsize(instate));
end
end

function [rppower_mean,rppower_peak,running_ave,pupil_ave] =...
    ripplepower_nostim(nss,nsstime,win,binsize,rippletime,...
    pupiltime,pupilsize,speedtime,speed,speedlimit)

rippletmp = rippletime(rippletime>=win(1) & rippletime<=win(2));
ripplehist = histcounts(rippletmp,win(1):binsize:win(2));

[~,~,bin] = histcounts(pupiltime,win(1):binsize:win(2));
pupilhist = cellfun(@(x) nanmean(pupilsize(bin==x)),num2cell(1:floor(diff(win)/binsize)));

[~,~,bin] = histcounts(speedtime,win(1):binsize:win(2));
speedhist = cellfun(@(x) nanmean(speed(bin==x)),num2cell(1:floor(diff(win)/binsize)));

[nssmean,nsspeak] = deal(cell(2,1));
for iRp = 1:2
    [~,~,bin] = histcounts(nsstime{iRp},win(1):binsize:win(2));
    nssmean{iRp} = cellfun(@(x) nanmean(nss{iRp}(bin==x)),num2cell(1:floor(diff(win)/binsize)));
    peaktemp = cellfun(@(x) max(nss{iRp}(bin==x)),num2cell(1:floor(diff(win)/binsize)),'UniformOutput',false);
    peaktemp(cellfun(@isempty,peaktemp)) = {NaN};
    nsspeak{iRp} = cell2mat(peaktemp);
end

pupilmobile = nanmean(pupilhist(abs(speedhist)>speedlimit));
[rppower_mean,rppower_peak] = deal(cell(4,1));
[running_ave,pupil_ave] = deal(NaN(4,1));
for iState = 1:4
    switch iState
        case 1 % ripple
            incondition = ripplehist>0;
        case {2,3} % immobile & no ripple & low-high level of pupil
            inconditiontemp = abs(speedhist)<speedlimit & ripplehist==0;
            if iState==2
                incondition = find(inconditiontemp & pupilhist<0.5*pupilmobile);
            elseif iState==3
                incondition = find(inconditiontemp & pupilhist>=0.5*pupilmobile);
            end
        case 4 % mobile
            incondition = find(abs(speedhist)>speedlimit);
    end
    
    for iRp = 1:length(nss)
        rppower_mean{iState} = cell2mat(cellfun(@(x) x(incondition),nssmean,'UniformOutput',false))';
        rppower_peak{iState} = cell2mat(cellfun(@(x) x(incondition),nsspeak,'UniformOutput',false))';
    end
    running_ave(iState) = mean(speedhist(incondition));
    pupil_ave(iState) = nanmean(pupilhist(incondition));
end
end