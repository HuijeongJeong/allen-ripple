clc; clearvars; close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

refList = {'medial';'lateral';'global'};

nclst = 3;
win_corr = 1;
% binsize = 0.01;
binsize = 0.02;
resolution = 1;

iCT = 1;
speedlimit = 2;

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

pair = [];
crosscorr = cell((nclst+2)*(nclst+1)/2,3,2);
clstidx = cell(nS,1);
ntime = NaN(nS,3);
for iS = 1:nS
    iS
    %% load data
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_ripples.mat'],'total_CA1_ripple','spontaneous_win');
    load([sdir(sessionList(iS)),'_Events.mat'],'running_speed','pupil_data');
    
    %% calculate index
    invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis));
    inthal = ismember(T.unit_id,tag.info.unit_id(tag.area.thalamus));
    
    [in,idx] = ismember(T.unit_id,unit_id.vis);
    clstidx{iS} = NaN(length(T.unit_id),1);
    clstidx{iS}(in) = cluster_idx.vis{nclst-1}(idx(in));
    clstidx{iS}(~in & invis) = 0;
    
    [~,idx] = ismember(T.unit_id,tag.info.unit_id);
    celltype = [tag.celltype.rs(idx), tag.celltype.fs(idx)];
    
    %% average psth for each cluster
    win = [spontaneous_win(1)-2*win_corr, spontaneous_win(2)+2*win_corr];
    spktime = win(1)+binsize/2:binsize:win(2)-binsize/2;
    nbin = length(spktime);
    
    spkhist = NaN(length(T.unit_id),nbin);
    spkhist(invis|inthal,:) =...
        cell2mat(cellfun(@(x) histcounts(x,win(1):binsize:win(2))*(1/binsize),...
        T.spike_time(invis|inthal),'UniformOutput',false));
    spkhistz = zscore(spkhist,[],2);
    
    data = NaN((nclst+1)*2+1,nbin);
    datact = NaN((nclst+1)*2+1,1);
    kClst = 1;
    for iClst = 1:nclst+2
        for iCT = 1:2
            if iClst==1
                inclst = inthal;
                if iCT==2
                    continue;
                end
            else
                inclst = invis & celltype(:,iCT) & clstidx{iS}==iClst-2;
                datact(kClst) = iCT;
            end
            data(kClst,:) = conv2(nanmean(spkhistz(inclst,:),1),...
                fspecial('Gaussian',[1 5*resolution],resolution),'same');
%             data(kClst,:) = nanmean(spkhistz(inclst,:),1);
            kClst = kClst+1;
        end
    end

    %% pupil 
    rippletime = cell2mat(total_CA1_ripple.ripple);
    rippletime = rippletime(:,1);
    rippletmp = rippletime(rippletime>=win(1) & rippletime<=win(2));
    ripplehist = histcounts(rippletmp,win(1):binsize:win(2))';
    list = find(ripplehist>0)';
    for i = list
        ripplehist(max([1,i-win_corr*1.5/binsize]):...
            min([length(ripplehist),i+win_corr*1.5/binsize])) = 1; % consider any bin <=1.5s apart from ripple as ripple window  
    end

    speedhist = interp1(running_speed.time,running_speed.velocity_conv,spktime');
    
    pupil_size = pupil_data.pupil_height/2.*pupil_data.pupil_width/2*pi;
    pupilhist = interp1(pupil_data.time,pupil_size,spktime');
    pupil_mobile = nanmean(pupilhist(abs(speedhist)>speedlimit));
    
    %% calculating correlation
    for iState = 1:3 %immobile & low ripple, immobile & high ripple, mobile
        if iState<3
            inanaltemp = abs(speedhist)<speedlimit & ripplehist==0;
            if iState==1
               inanal = inanaltemp & pupilhist<pupil_mobile*0.5; 
            else
               inanal = inanaltemp & pupilhist>=pupil_mobile*0.5;
            end
        else
            inanal = abs(speedhist)>=speedlimit;
        end
        
        analwin = [find(diff([0;inanal])==1),find(diff([inanal;0])==-1)];
        analwin(diff(analwin,[],2)<3/binsize,:) = []; % use only window that is longer than 3s
        ntime(iS,iState) = sum(diff(analwin,[],2)*binsize); 
        if ntime(iS,iState)<50 % if total time of any state is below than 50s, exclude the session from analysis
            continue;
        end
        
        for iPair = 1 % rs-rs, fs-fs
            datasub = data(datact==iPair | isnan(datact),:);
            
            k = 1;
            for i = 1:nclst+1
                for j = i+1:nclst+2
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
                    crosscorr{k,iState,iPair}(iS,:) = nanmean(corrtemp,1);
                    k = k+1;
                end
            end
        end
    end
end


%%
nmod = cell2mat(cellfun(@(x) [sum(x==1),sum(x==2),sum(x==3)],clstidx,'UniformOutput',false));
out = sum(nmod==0,2)>0;
in = ~out & sum(ntime>50,2)==3;

%%
iState = 1;
iCT = 1;

timecorr = -win_corr:binsize:win_corr;
ct = [cbrewer('qual','Set2',3);cbrewer('qual','Dark2',3)];
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20/6 20/6]);
hold on;
for iClst = 1:3
    plot(timecorr,crosscorr{iClst+1,iState,iCT},...
        'Color',ct(iClst,:),'LineWidth',0.35);
end
plot(timecorr,nanmean(crosscorr{1,iState,iCT}),...
       'Color','k','LineWidth',1);
for iClst = 1:3
   plot(timecorr,nanmean(crosscorr{iClst+1,iState,iCT}),...
       'Color',ct(iClst+3,:),'LineWidth',1); 
end
plot([0 0],[-0.4 0.6],'k:')
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,...
    'XTick',-1:0.5:1,'YTick',-0.4:0.2:0.6,'YLim',[-0.4 0.6]);
xlabel('Lag (s)');
ylabel('Correlation');
cd('D:\OneDrive - University of California, San Francisco\figures\2.allen\fig4');
print(fHandle,'-depsc','-painters','crosscorr_thal.ai')

%%
inpre = timecorr>=-0.5 & timecorr<0;
inpost = timecorr>0 & timecorr<=0.5;

ct = [0.8 0.8 0.8;cbrewer('qual','Set2',3);0 0 0;cbrewer('qual','Dark2',3)];

selec = NaN(nS,4);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3.5 20*0.18]);
j = 0;
hold on;
for iClst = [0,1,2,3]
    auc = NaN(nS,2);
    for iS = 1:nS
         auc(iS,1) = trapz(timecorr(inpre),crosscorr{iClst+1,iState,iCT}(iS,inpre));
        auc(iS,2) = trapz(timecorr(inpost),crosscorr{iClst+1,iState,iCT}(iS,inpost));
%         auc(iS,1) = trapz(timecorr(inpre),crosscorr{iClst+1,iState,iCT}(iS,inpre));
%         auc(iS,2) = trapz(timecorr(inpost),crosscorr{iClst+1,iState,iCT}(iS,inpost));
    end
    [~,ptemp,~,stat] = ttest(abs(auc(:,2)),abs(auc(:,1)));
    p(iClst+1) = ptemp*4;
    t(iClst+1) = stat.tstat;
    selec(:,iClst+1) = abs(auc(:,1))-abs(auc(:,2));
    scatter(rand(size(auc,1),1)*0.5+j,selec(:,iClst+1),3,ct(iClst+1,:),'filled');
    errorbar(0.25+j,nanmean(selec(:,iClst+1)),...
        nanstd(selec(:,iClst+1))/sqrt(sum(~isnan(selec(:,iClst+1)))),'Color',ct(iClst+5,:),'LineWidth',1);
    j = j+1;
end
plot([-0.5 4],[0 0],'k:')
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',0.25:1:3.25,...
    'YTick',-0.12:0.06:0.06,'YLim',[-0.13 0.09],'XLim',[-0.5 4],'XTickLabel',...
    {'Nomod';'iAct';'dAct';'Inh'},'XTickLabelRotation',45);
ylabel('\Deltaabsolute AUC (Pre-Post)');

cd('D:\OneDrive - University of California, San Francisco\figures\2.allen\fig4');
print(fHandle,'-depsc','-painters','crosscorr_thal_auc.ai')

%%
[p,t] = deal(nan(3,1));
for iClst = 1:3
    [~,ptemp,~,stat]=ttest(selec(~isnan(selec(:,1)),1),selec(~isnan(selec(:,1)),iClst+1));
    p(iClst) = ptemp*3;
    t(iClst) = stat.tstat;
end

% [tbl,rm] = simple_mixed_anova(selec)

% %%
% fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2 20*0.18]);
% % j = 0;
% hold on;
% for iClst = [1,3]
%     auc = NaN(nS,2);
%     for iS = 1:nS
%          auc(iS,1) = trapz(timecorr(inpre),crosscorr{iClst+1,iState,iCT}(iS,inpre));
%         auc(iS,2) = trapz(timecorr(inpost),crosscorr{iClst+1,iState,iCT}(iS,inpost));
% %         auc(iS,1) = trapz(timecorr(inpre),crosscorr{iClst+1,iState,iCT}(iS,inpre));
% %         auc(iS,2) = trapz(timecorr(inpost),crosscorr{iClst+1,iState,iCT}(iS,inpost));
%     end
%     plot(1:2,abs(auc),'Color',ct(iClst,:),'LineWidth',0.35);
% errorbar(1:2,nanmean(abs(auc)),nanstd(abs(auc))/sqrt(sum(~isnan(auc(:,1)))),...
%     'Color',ct(iClst+3,:),'LineWidth',1,'CapSize',3);
% %     [~,p(iClst)] = ttest(auc(:,2),auc(:,1));
% %     selec = auc(:,1)./auc(:,2);
% %     scatter(rand(size(auc,1),1)*0.5+j,selec,3,ct(iClst,:),'filled');
% %     errorbar(0.25+j,nanmean(selec),...
% %         nanstd(selec)/sqrt(sum(~isnan(selec))),'Color',ct(iClst+3,:),'LineWidth',1);
% %     j = j+1;
% end
% % plot([0.5 2.5],[0 0],'k:')
% set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',[1 2],...
%     'YTick',-0.2:0.1:0.2,'YLim',[0 0.2],'XLim',[0.5 2.5],'XTickLabel',{'Pre';'Post'},...
%     'XTickLabelRotation',45);
% ylabel('AUC');
% 
% 
% cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig4');
% print(fHandle,'-depsc','-painters','crosscorr_thal_auc.ai')
%%
iClst = 1;

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20/6 20/6]);
hold on;
for iClst = 1:3
%     
plot(timecorr,flip(crosscorr{pair(:,1)==iClst+1&pair(:,2)==5,iState,iCT},2),...
    'Color',ct(iClst,:),'LineWidth',0.35);

plot(timecorr,nanmean(flip(crosscorr{pair(:,1)==iClst+1&pair(:,2)==5,iState,iCT},2)),...
    'Color',ct(iClst+4,:),'LineWidth',1);
end
plot([0 0],[-0.4 0.6],'k:')
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,...
    'XTick',-1:0.5:1,'YTick',-0.4:0.2:0.6,'YLim',[-0.4 0.6]);
xlabel('Lag (s)');
ylabel('Correlation');
cd('D:\OneDrive - University of California, San Francisco\figures\2.allen\fig4');
print(fHandle,'-depsc','-painters','crosscorr_inh.ai')

%%
[p,t]=deal(nan(3,1));
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 20*0.18]);
hold on;
selec = NaN(nS,3);
for iClst = 1:3
data = flip(crosscorr{pair(:,1)==iClst+1&pair(:,2)==5,iState,iCT},2);
auc = NaN(nS,2);
for iS = 1:nS
    auc(iS,1) = trapz(timecorr(inpre),data(iS,inpre));
    auc(iS,2) = trapz(timecorr(inpost),data(iS,inpost));
end
[~,ptemp,~,stat] = ttest(abs(auc(:,2)),abs(auc(:,1)));
    p(iClst+1) = ptemp*3;
    t(iClst+1) = stat.tstat;
selec(:,iClst) = abs(auc(:,1))-abs(auc(:,2));
    scatter(rand(size(auc,1),1)*0.5+iClst-1,selec(:,iClst),3,ct(iClst,:),'filled');
    errorbar(0.25+iClst-1,nanmean(selec(:,iClst)),...
        nanstd(selec(:,iClst))/sqrt(sum(~isnan(selec(:,iClst)))),'Color',ct(iClst+4,:),'LineWidth',1);
end
    plot([-0.5 3],[0 0],'k:')
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',0.25:1:2.25,...
        'YTick',-0.12:0.06:0.06,'YLim',[-0.13 0.09],'XLim',[-0.5 3],'XTickLabel',{'Nomod';'iAct';'dAct'},...
    'XTickLabelRotation',45);
ylabel('\Deltaabsolute AUC (Pre-Post)');

cd('D:\OneDrive - University of California, San Francisco\figures\2.allen\fig4');
print(fHandle,'-depsc','-painters','crosscorr_inh_auc.ai')

%%
[p,t] = deal(nan(2,1));
for iClst = 1:2
    [~,ptemp,~,stat]=ttest(selec(~isnan(selec(:,1)),1),selec(~isnan(selec(:,1)),iClst+1));
    p(iClst) = ptemp*2;
    t(iClst) = stat.tstat;
end