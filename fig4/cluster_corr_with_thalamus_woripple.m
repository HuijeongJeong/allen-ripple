clearvars; clc; close all;

load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

binsize = 2;
speedlimit = 2;
threshold = [2 5];   %threshold for ripple beginning/end and peak
duration = [30 20 250];    %min inter-ripple inter, min ripple length, max ripple length
Fs = 1250;

typeList = {'rs';'fs'};
celltype = [tag.celltype.rs tag.celltype.fs];
iT = 1;

cidx_vis = deal(cell(nS,1));
[ncorr_sp,ncorr_dg,scorr_dg] = deal(cell(nS,4));
[structure_vis,structure_thal] = deal(cell(nS,1));
[pupil_sp,speed_sp,nbin_sp] = deal(NaN(nS,4));

%%
for iS = 1:nS
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_Events.mat'],'running_speed','pupil_data','drifting_gratings');
    load([sdir(sessionList(iS)),'_ripples.mat'],'filtered_lfp','total_CA1_ripple','spontaneous_win');
    load([sdir(sessionList(iS)),'_stimulus_evoked.mat'],'fr_gratings','beh_gratings');
    
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
    dglist = sort(unique(drifting_gratings.stimulus_condition_id));
    ndg = length(dglist);
    dgIndex = false(length(drifting_gratings.stimulus_condition_id),ndg);
    for idg = 1:ndg
        dgIndex(:,idg) = drifting_gratings.stimulus_condition_id==dglist(idg);
    end
    %%
    invis = ismember(T.unit_id,tag.info.unit_id(celltype(:,iT) & tag.area.vis));
    inthal = ismember(T.unit_id,tag.info.unit_id(tag.area.thalamus));
    
    unitid_tmp = [T.unit_id(invis); T.unit_id(inthal)];
    areaidx = false(sum(invis|inthal),2);
    areaidx(1:sum(invis),1) = true;
    areaidx(sum(invis)+1:end,2) = true;
        
    nCell = sum(areaidx);
    if sum(nCell==0)>0
        continue;
    end
    
    k = 3;
    [in,idxx] = ismember(unitid_tmp(areaidx(:,1)),unit_id.vis);
    cidx_vis{iS} = zeros(length(idxx),1);
    cidx_vis{iS}(in) = cluster_idx.vis{k-1}(idxx(in));
    
    [~,idxx] = ismember(unitid_tmp(areaidx(:,1)),T.unit_id);
    structure_vis{iS} = T.ecephys_structure_acronym(idxx);
    [~,idxx] = ismember(unitid_tmp(areaidx(:,2)),T.unit_id);
    structure_thal{iS} = T.ecephys_structure_acronym(idxx);
    
    pupil_size = pupil_data.pupil_height/2.*pupil_data.pupil_width/2*pi;
  
    ripples = cell2mat(ripples);
    dg_ripple = cellfun(@(x) sum(ripples(:,1)>=x(1) & ripples(:,1)<=x(2)),...
        mat2cell(drifting_gratings.window,ones(size(drifting_gratings.window,1),1),2));
    for iState = 1:4
        fprintf('%sth block of %sth session\n',num2str(iState),num2str(iS));
        [~,~,ncorr_sp{iS,iState},pupil_sp(iS,iState),speed_sp(iS,iState),nbin_sp(iS,iState)] =...
            noisecorr_nostim([T.spike_time(invis);T.spike_time(inthal)],...
            spontaneous_win,binsize,unitid_tmp,areaidx,ripples(:,1),iState,...
            pupil_data.time,pupil_size,running_speed.time,running_speed.velocity_conv,speedlimit);
        
        [~,~,scorr_dg{iS,iState},ncorr_dg{iS,iState},~,~] =...
            paircorr([fr_gratings.spikeTime(invis);fr_gratings.spikeTime(inthal)],...
            [0 2],binsize,unitid_tmp,areaidx,dg_ripple,iState,dgIndex,...
            cell2mat(beh_gratings.velocity_ave'),speedlimit,cell2mat(beh_gratings.pupil_area_ave'));


    end
end

%%
c = cbrewer('div','RdBu',50);
c([1:5, 46:50],:) = [];
thalList = {'LG';'LP';'MG'};

stateList = {'Immobile-low pupil';'Immobile-high pupil';'Mobile'};

nmod = cell2mat(cellfun(@(x) [sum(x==1),sum(x==2),sum(x==3)],cidx_vis,'UniformOutput',false));
no_ripple_mod = sum(nmod==0,2)>0;
in = sum(nbin_sp(:,2:4)>50/binsize,2)==3 & ~no_ripple_mod;
%%

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20/2 20/6]);
for iState = 2:4
   axes('Position',axpt(3,10,iState-1,1:9));
   m_sp = NaN(nS,4,4);
   for iR = 1:4
      switch iR
          case {1,2,3}
              inthal = cellfun(@(x) contains(x,thalList{iR}),structure_thal(in),'UniformOutput',false);
          case 4
              inthal = cellfun(@(x) ~contains(x,thalList),structure_thal(in),'UniformOutput',false);
      end
      no_thal = cellfun(@sum,inthal)<5;
      for iClst = 1:k+1
          m_sp(in==1,iR,iClst) = cellfun(@(x,y,z) nanmean(nanmean(x(z==iClst-1,y),2)),...
              ncorr_sp(in==1,iState),inthal,cidx_vis(in==1));
      end
      m_sp(no_thal,iR,:) = NaN;
   end
   imagesc(squeeze(nanmean(m_sp,1)));
   set(gca,'CLim',[-0.1 0.1],'Box','off','TickDir','out','FontSize',7,...
       'XTick',1:4,'XTickLabel',{'Nomod','iAct','dAct','Inh'},'XTickLabelRotation',45,...
       'YTick',1:4,'YTickLabel',{'LGN','LP','MGN','Other'});
   if iState>2
       set(gca,'YTickLabel',[]);
   end
   title(stateList{iState-1});
end
colormap(c);
% cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig4');
% print(fHandle,'-depsc','-painters',...
%     ['correlation_',typeList{iT},'_',num2str(binsize*1000),'bin_wo_ripple_thal.ai'])

%%
ct = [cbrewer('qual','Dark2',3);cbrewer('qual','Set2',3)];
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.15 20*0.18]);
hold on;
nc_sp = NaN(nS,3,3);
for iClst = 1:k
    for iState = 2:4
        nc_sp(in,iState-1,iClst) = cellfun(@(x,y) nanmean(nanmean(x(y==iClst,:),2)),...
            ncorr_sp(in,iState),cidx_vis(in));
    end
    plot(1:3,squeeze(nc_sp(:,:,iClst)),'Color',ct(iClst+3,:),'LineWidth',0.35);
end
for iClst = 1:k
    errorbar(1:3,nanmean(squeeze(nc_sp(:,:,iClst))),nanstd(squeeze(nc_sp(:,:,iClst)))/...
        sqrt(sum(in)),'Color',ct(iClst,:),'CapSize',3);
end
plot([0 4],[0 0],'k:');
xlim([0.5 3.5]);
ylim([-0.15 0.2]);
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'YTick',-0.1:0.1:0.2,...
    'XTick',1:3,'XTickLabel',{'low pupil';'high pupil';'mobile'},'XTickLabelRotation',45);
ylabel('Noise correlation');

% cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig4');
% print(fHandle,'-depsc','-painters','noise_corr_thal_woripple.ai');

tbl = simple_mixed_anova(nc_sp);


%% comparison b/w v1-lgn vs. v1-lp
thalList2 = {'LGd';'LP';'MG'};

nc_sp = NaN(nS,3,4,4);
for iVis = 1:2
    if iVis==1
        invis = cellfun(@(x) strcmp(x,'VISp'),structure_vis(in),'UniformOutput',false);
    else
        invis = cellfun(@(x) ~strcmp(x,'VISp'),structure_vis(in),'UniformOutput',false);
    end
    for iR = 1:2
        inthal = cellfun(@(x) contains(x,thalList2{iR}),structure_thal(in),'UniformOutput',false);
        no_thal = cellfun(@sum,inthal)<3;
        for iClst = 1:k+1
            inclst = cellfun(@(x) x==iClst-1,cidx_vis(in),'UniformOutput',false);
            no_vis = cellfun(@(x,y) sum(x&y),inclst,invis)<3;
            for iState = 2:4
                nc_temp = cellfun(@(x,y,z,w) nanmean(x(y&w,z),2),ncorr_sp(in,iState),...
                    invis,inthal,inclst,'UniformOutput',false);
                nc_temp(cellfun(@isempty,nc_temp)) = {nan};
                nc_sp(in,iState-1,iClst,iR+(iVis-1)*2) = cellfun(@nanmean,nc_temp);
            end
        end
    end
end

dir = 'D:\OneDrive - University of California, San Francisco\figures\allen';
vislist = {'V1';'HVA'};
clusterlist = {'iAct';'dAct';'Inh'};
clr = [0,0,0;cbrewer('qual','Dark2',3)];
p = nan(3,2,2);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 6]);
for iState = 1
    for iClst = 2:4
        for i =1:2
            axes('Position',axpt(3,2,iClst-1,i,axpt(10,10,3:10,1:9))); hold on;
            in = sum(squeeze(isnan(nc_sp(:,iState,iClst,[1:2]+(i-1)*2))),2)==0;
            sum(in)
            data = squeeze(nc_sp(in,iState,iClst,[1:2]+(i-1)*2));
            plot(1:2,data,'Color',[0.6 0.6 0.6],'LineWidth',0.35);
            plot(1:2,mean(data),'k','LineWidth',0.75);
            errorbar(i,mean(data(:,i)),std(data(:,i))/sqrt(sum(in)),'r','CapSize',3);
            errorbar(2/i,mean(data(:,2/i)),std(data(:,2/i))/sqrt(sum(in)),'k','CapSize',3);
            plot([0 5],[0 0],'k:','LineWidth',0.35);
            ylim([-0.25 0.4])
            xlim([0.5 2.5]);
            [~,p,~,stat] = ttest(data(:,1),data(:,2));
            ppair(iClst,i) = p;
            statpair(iClst,i) = stat.tstat;
            
            for ii = 1:2
               [~,p,~,stat] = ttest(data(:,ii));
               psig(iClst,i,ii) = p*2;
               statsig(iClst,i,ii) = stat.tstat;
            end
%             tbl = simple_mixed_anova(data);
%             p(iClst-1,i,:) = tbl.pValue([1,3]);
%             if ppair(iClst,i)<0.05
%                text(1.5,0.35,'*','FontSize',8); 
%             end
%             if p(iClst-1,i,1)<0.05
%                text(2.3,0,'*','FontSize',8); 
%             end
            set(gca,'Box','off','TickDir','out','FontSize',8,'XTick',1:2,...
                'YTick',-0.2:0.2:0.4);
            if i==1
                set(gca,'XTickLabel',[]);
                title(clusterlist{iClst-1});
            else
                set(gca,'XTickLabel',{'LGN','LP'},'XTickLabelRotation',45);
            end
            if iClst==2
               ylabel({'Spont. correlation';['(',vislist{i},' vs. thalamus)']});
            else
                set(gca,'YTickLabel',[]);
            end
        end
    end
end
dir = 'D:\OneDrive - University of California, San Francisco\figures\2.allen\figS7';
print(fHandle,'-depsc','-painters',[dir,'\noise_corr_thal_woripple_subvis.ai']);

tbl = simple_mixed_anova(squeeze(nc_sp(in,iState,iClst,[1:2]+(i-1)*2)));

%%

function [pairid_1,pairid_2,signalcorr,noisecorr,ntrial,pupil] =...
    paircorr(spikeTime,win,binsize,unit_id,pairIdx,trialripple,state,...
    condition,speed,speedlimit,pupilsize)

spkhist = cellfun(@(x) cell2mat(cellfun(@(y) histcounts(y,win(1):binsize:win(2)),x,...
    'UniformOutput',false)),spikeTime,'UniformOutput',false);

nCell = length(spikeTime);
in = false(nCell,nCell);
in(pairIdx(:,1),pairIdx(:,2)) = true;
in = in(:);

x = repmat(unit_id,1,length(unit_id));
y = repmat(unit_id',length(unit_id),1);

pairid = [x(:),y(:)];
pairid = pairid(in,:);

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
pairid_1 = reshape(pairid(:,1),sum(pairIdx(:,1)),sum(pairIdx(:,2)));
pairid_2 = reshape(pairid(:,2),sum(pairIdx(:,1)),sum(pairIdx(:,2)));

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
    signalcorr = scorr(in);
    signalcorr = reshape(signalcorr,sum(pairIdx(:,1)),sum(pairIdx(:,2)));
    
    noise = cell2mat(noise);
    ncorr = corr(noise');
    ncorr = ncorr(:);
    noisecorr = ncorr(in);
    noisecorr = reshape(noisecorr,sum(pairIdx(:,1)),sum(pairIdx(:,2)));
    
    ntrial = sum(instate);
    pupil = nanmean(pupilsize(instate));
else
    ntrial = 0;
    pupil = NaN;
    [signalcorr,noisecorr] = deal(NaN(sum(pairIdx(:,1)),sum(pairIdx(:,2))));
end

end


function [pairid_1,pairid_2,noisecorr,avepupil,avespeed,nbin] =...
    noisecorr_nostim(spikeTime,win,binsize,unit_id,pairIdx,rippletime,...
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
    case {2,3} % immobile & no ripple & low-middle-high level of pupil
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

nCell = length(spikeTime);
in = false(nCell,nCell);
in(pairIdx(:,1),pairIdx(:,2)) = true;
in = in(:);

x = repmat(unit_id,1,length(unit_id));
y = repmat(unit_id',length(unit_id),1);

pairid = [x(:),y(:)];
pairid = pairid(in,:);

if ~isempty(incondition)
    spkhist = spkhist(:,incondition);
    noisecorr = corr(spkhist');
    noisecorr = noisecorr(:);
    noisecorr = noisecorr(in);
    
    pairid_1 = reshape(pairid(:,1),sum(pairIdx(:,1)),sum(pairIdx(:,2)));
    pairid_2 = reshape(pairid(:,2),sum(pairIdx(:,1)),sum(pairIdx(:,2)));
    noisecorr = reshape(noisecorr,sum(pairIdx(:,1)),sum(pairIdx(:,2)));
    
    nbin = length(incondition);
    
    avepupil = nanmean(pupilhist(incondition));
    avespeed = nanmean(speedhist(incondition));
else
    nbin = 0;
    [avepupil,avespeed] = deal(NaN);
    noisecorr = NaN(sum(pairIdx(:,1)),sum(pairIdx(:,2)));
end
end