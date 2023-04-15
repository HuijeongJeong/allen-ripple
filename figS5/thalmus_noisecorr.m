clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

arealist = {'LG';'LP';'MG'};
nA = length(arealist);

noisecorr_sp = cell(nS,nA*(nA+1)/2+nA+1);
[time_sp,pupil_sp,nbin_sp] = deal(NaN(nS,4));
[areaidx,distance,withinarea] = deal(cell(nS,1));

speedlimit = 2;
binsize = 2;

threshold = [2 5];   %threshold for ripple beginning/end and peak
duration = [30 20 250];    %min inter-ripple inter, min ripple length, max ripple length
Fs = 1250;

%%
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

    inthal = ismember(T.unit_id,tag.info.unit_id(tag.area.thalamus));
    nCell = sum(inthal);
    
    if nCell==0 
        continue;
    end
    
    unitid_tmp = T.unit_id(inthal);
    [~,idx] = ismember(unitid_tmp,tag.info.unit_id);
    structure = tag.info.structure(idx);
    ccf = tag.info.ccf(idx,:);
    area = tag.info.structure(idx);
    
    areaidx{iS} = zeros(nCell,1);
    for iA = 1:nA
        areaidx{iS}(contains(area,arealist{iA})) = iA; 
    end

    %%    
    
    pupil_size = pupil_data.pupil_height/2.*pupil_data.pupil_width/2*pi;

    noise_corr_sp = NaN(nCell*(nCell-1)/2,4);
    ripples = cell2mat(ripples);
    
    for iState = 1:4
        [pairid,noise_corr_sp(:,iState),time_sp(iS,iState),pupil_sp(iS,iState),nbin_sp(iS,iState)] =...
            noisecorr_nostim(T.spike_time(inthal),spontaneous_win,binsize,unitid_tmp,ripples(:,1),...
            iState,pupil_data.time,pupil_size,running_speed.time,running_speed.velocity_conv,speedlimit);
    end
    ccfout = sum(ccf<0,2)>0;
    ccf(ccfout,:) = NaN;
    pair_distance = pdist(ccf)';

    [~,idx] = ismember(pairid,unitid_tmp);
    aidxtmp = areaidx{iS}(idx);
    areatmp = structure(idx);
    kClst = 1;
    for iA = 1:nA+1
        for jA = iA:nA+1
            in = aidxtmp(:,1)==iA-1 & aidxtmp(:,2)==jA-1;

            noisecorr_sp{iS,kClst} = noise_corr_sp(in,:);
            distance{iS,kClst} = pair_distance(in);
            withinarea{iS,kClst} = cellfun(@(x,y) strcmp(x,y),areatmp(in,1),areatmp(in,2));
            kClst = kClst+1;
        end
    end
end

%%
nmod = cell2mat(cellfun(@(x) [sum(x==0),sum(x==1),sum(x==2),sum(x==3)],areaidx,'UniformOutput',false));
out = sum(nmod==0,2)>0;

x = [1 1 1 1 2 2 2 3 3 4];
y = [1 2 3 4 2 3 4 3 4 4];

stateList = {'all';'Immobile-low pupil';'Immobile-high pupil';'Mobile'};
c = cbrewer('div','RdBu',50);
c([1:5, 46:50],:) = [];

%%
data = NaN(nS,length(x),4);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20/2 20/6]);
for iState = 2:4
    in = sum(nbin_sp(:,2:4)>50/binsize,2)==3;
    data(in,:,iState) = cellfun(@(x) nanmean(x(:,iState)),noisecorr_sp(in,:));
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
        'XTickLabel',{'LGN';'LP';'MGN';'Other'},'YTick',1:4,'YTickLabel',...
        {'LGN';'LP';'MGN';'Other'},'XTickLabelRotation',45,'Box','off','TickDir','out');
    title(stateList{iState});
    if iState>2
        set(gca,'YTickLabel',[]);
    end
end
% cd('D:\OneDrive - University of California, San Francisco\figures\allen\figS5');
% print(fHandle,'-depsc','-painters',...
%     ['correlation_',num2str(binsize*1000),'bin_wo_ripple_thalamus.ai'])


%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.19 20/6]);
iState = 1;
in = sum(nbin_sp(:,2:4)>50/binsize,2)==3 & ~out;
data(in,:,iState) = cellfun(@(x) nanmean(x(:,iState)),noisecorr_sp(in,:));

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
    'XTickLabel',{'LGN';'LP';'MGN';'Other'},'YTick',1:4,'YTickLabel',...
        {'LGN';'LP';'MGN';'Other'},'XTickLabelRotation',45,'Box','off','TickDir','out');
% print(fHandle,'-depsc','-painters',...
%     ['correlation_',num2str(binsize*1000),'bin_wo_ripple_thalamus_all.ai'])


%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.15 4]);
hold on;
in = sum(isnan(squeeze(nanmean(data(:,x~=y,2:4),2))),2)==0;
mwithin = squeeze(nanmean(data(in,x==y,2:4),2));
mbetween = squeeze(nanmean(data(in,x~=y,2:4),2));
plot(1:3,mbetween,'Color',[0.8 0.8 0.8]);
plot(1:3,mwithin,'Color',[1 0.8 0.8]);
errorbar(1:3,nanmean(mbetween),nanstd(mbetween)/sqrt(size(mbetween,1)),'k')
errorbar(1:3,nanmean(mwithin),nanstd(mwithin)/sqrt(size(mwithin,1)),'r')
plot([0.5 3.5],[0 0],'k:','LineWidth',0.35);
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',1:3,...
    'XTickLabel',{'Low pupil';'High pupil';'Mobile'},...
    'XTickLabelRotation',45,'YTick',-0.1:0.1:0.2)
ylim([-0.1 0.25])
xlim([0.5 3.5])
ylabel('Correlation')
tbl = simple_mixed_anova(cat(3,mbetween,mwithin));
p = tbl.pValue([5,3,7]); %state, area, state*area
for i = 1:3
   text(0.7,0.24-0.05*(i-1),['P = ',num2str(p(i))],'FontSize',6)
end
% print(fHandle,'-depsc','-painters',...
%     ['correlation_',num2str(binsize*1000),'bin_wo_ripple_thalamus_w_state.ai'])
%%
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