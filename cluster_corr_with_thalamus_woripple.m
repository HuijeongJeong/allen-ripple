clearvars; clc; close all;

load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

binsize = 2;
speedlimit = 2;

typeList = {'rs';'fs'};
celltype = [tag.celltype.rs tag.celltype.fs];
iT = 1;

cidx_vis = deal(cell(nS,1));
ncorr_sp = deal(cell(nS,3));
structure_vis = deal(cell(nS,1));
[ntrial_dg,time_sp] = deal(NaN(nS,3,2));

for iS = 1:nS
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_Events.mat'],'running_speed');
    load([sdir(sessionList(iS)),'_ripples.mat'],'total_CA1_ripple','spontaneous_anal_win','spontaneous_win');
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
    
    win_sp = cell(3,1);
    spontaneous_mobile = [spontaneous_win(1), spontaneous_anal_win(1,1);...
        spontaneous_anal_win(1:end-1,2), spontaneous_anal_win(2:end,1);...
        spontaneous_anal_win(end,1) spontaneous_win(2)];
    spontaneous_mobile(diff(spontaneous_mobile,[],2)<=binsize,:) = [];
    win_sp{1} = spontaneous_anal_win;
    win_sp{2} = spontaneous_anal_win;
    win_sp{3} = spontaneous_mobile;
    
    ripples = sortrows(cell2mat(total_CA1_ripple.ripple));
    for iState = 1:3
        fprintf('%sth block of %sth session\n',num2str(iState),num2str(iS));
        [~,~,ncorr_sp{iS,iState},time_sp(iS,iState)] =...
            noisecorr_nostim([T.spike_time(invis);T.spike_time(inthal)],...
            win_sp{iState},binsize,unitid_tmp,areaidx,ripples(:,1),iState);
    end
end


nmod = cell2mat(cellfun(@(x) [sum(x==1),sum(x==2),sum(x==3)],cidx_vis,'UniformOutput',false));
out = sum(nmod==0,2)>0;

clr = {[0.6 0.6 0.6],[0 0 0]; [0.6 0.6 1],[0 0 1]; [1 0.6 0.6],[1 0 0]};

iCT = 1;
iState = 1;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 3]);
    hold on;
    nc_sp = NaN(nS,k+1);
    for iClst = 1:k+1
        nc_sp(:,iClst) = cellfun(@(x,y,z) nanmean(nanmean(x(y==iClst-1,:),2)),...
            ncorr_sp(:,iState),cidx_vis);
    end
    plot(1:4,nc_sp,'Color',clr{1,1});
    errorbar(1:4,nanmean(nc_sp),nanstd(nc_sp)/sqrt(size(nc_sp,1)),'k');
    plot([0 5],[0 0],'k:');
    xlim([0 5]);
    ylim([-0.12 0.12]);
    set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,'YTick',-0.1:0.05:0.1,...
        'XTick',1:4,'XTickLabel',{'nomod';'Iact';'Dact';'Inh'},'XTickLabelRotation',45);
        ylabel('Noise correlation');

cd('D:\OneDrive - University of California, San Francisco\figures\allen\correlation')
print(fHandle,'-dtiff','-r600','noise_corr_thal_woripple.tif');



iHp = 1;
ylimit = [-0.1 0.2; -0.2 0.3];
for iRef = 1:2
for iHT = 1:2
    for iLoc = 1:3
        [nc_dg,sc_dg,nc_sp] = deal(NaN(nS,3,k+1));
        
        inthal = cellfun(@(x,y,z) strcmp(x,hippoList{iHp}) & y(:,iHT) & z(:,iLoc),...
            structure_hippo,celltype_hippo,locIdx,'UniformOutput',false);
%         inhippo = cellfun(@(x,y,z) strcmp(x,hippoList{iHp}) & y(:,iHT) & z(:,iLoc),...
%             structure_hippo(~out),celltype_hippo(~out),locIdx(~out),'UniformOutput',false);
        for iClst = 1:k+1
            nc_dg(:,:,iClst) = cellfun(@(x,y,z) nanmean(nanmean(x(y==iClst-1,z),2)),...
                ncorr_dg(:,:,iRef),repmat(cidx_vis,1,3),repmat(inthal,1,3));
            nc_sp(:,:,iClst) = cellfun(@(x,y,z) nanmean(nanmean(x(y==iClst-1,z),2)),...
                ncorr_sp(:,:,iRef),repmat(cidx_vis,1,3),repmat(inthal,1,3));
            sc_dg(:,:,iClst) = cellfun(@(x,y,z) nanmean(nanmean(x(y==iClst-1,z),2)),...
                scorr_dg(:,:,iRef),repmat(cidx_vis,1,3),repmat(inthal,1,3));
        end
        
        fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 5]);
        for iD = 1:2
            switch iD
                case 1
                    s = sc_dg(~out,:,:); n = nc_dg(~out,:,:);
                case 2
                    n = nc_sp(~out,:,:);
            end
            
            if iD==1
            axes('Position',axpt(2,2,iD,1,axpt(10,10,2:10,1:9)))
            hold on;
            for iState = 1:3
                plot(1:k+1,squeeze(s(:,iState,:)),'Color',clr{iState,1},'LineWidth',0.3);
            end
            for iState = 1:3
                errorbar(1:k+1,nanmean(squeeze(s(:,iState,:))),nanstd(squeeze(s(:,iState,:)))./...
                    sqrt(sum(~isnan(squeeze(s(:,iState,:))))),'Color',clr{iState,2},'LineWidth',0.5,'CapSize',3);
            end
            plot([0 5],[0 0],'k:','LineWidth',0.35);
            ylim(ylimit(iHT,:))
            xlim([0 5])
            set(gca,'XTick',1:4,'XTickLabel',[],'Box','off','TickDir','out','FontSize',4,...
                'LineWidth',0.35,'YTick',ylimit(iHT,1):abs(ylimit(iHT,1)):ylimit(iHT,2));
            ylabel('Signal correlation','FontSize',5);
            title('Gratings','FontSize',5);
            end
            
            axes('Position',axpt(2,2,iD,2,axpt(10,10,2:10,1:9)))
            hold on;
            for iState = 1:3
                plot(1:k+1,squeeze(n(:,iState,:)),'Color',clr{iState,1},'LineWidth',0.3);
            end
            for iState = 1:3
                errorbar(1:k+1,nanmean(squeeze(n(:,iState,:))),nanstd(squeeze(n(:,iState,:)))./...
                    sqrt(sum(~isnan(squeeze(n(:,iState,:))))),'Color',clr{iState,2},'LineWidth',0.5,'CapSize',3);
            end
            plot([0 5],[0 0],'k:','LineWidth',0.35);
            ylim(ylimit(iHT,:))
            xlim([0 5])
            set(gca,'XTick',1:4,'XTickLabel',{'Nomod';'I act';'D act';'Inh'},'XTickLabelRotation',45,...
                'YTick',ylimit(iHT,1):abs(ylimit(iHT,1)):ylimit(iHT,2),'Box','off','TickDir',...
                'out','FontSize',4,'LineWidth',0.35);
            if iD==1
                ylabel('Noise correlation','FontSize',5);
            else
                set(gca,'YTickLabel',[]);
                    title('Spontaneous','FontSize',5);
            end
        end
        print(fHandle,'-dtiff','-r600',...
            ['D:\heejeong\OneDrive\Research\DataFig\1.allen-andermann\clustering\corr\pair_corr_with_',...
            hippoList{iHp},typeList{iHT},'_',typeList{iT},locList{iLoc},'_wo_',refList2{iRef},'ripple.tif']);
        close all
    end
end
end


function [pairid_1,pairid_2,noisecorr,time] = noisecorr_nostim(spikeTime,win,binsize,unit_id,pairIdx,rippletime,iRipple)

nWin = size(win,1);
[spkhist,ripplehist] = deal(cell(1,nWin));
for iW = 1:nWin
    if diff(win(iW,:))<binsize
        continue;
    end
    spiketmp = cellfun(@(x) x(x>=win(iW,1) & x<=win(iW,2)),spikeTime,'UniformOutput',false);
    spkhist{iW} = cell2mat(cellfun(@(x) histcounts(x,win(iW,1):binsize:win(iW,2)),spiketmp,'UniformOutput',false));
    rippletmp = rippletime(rippletime>=win(iW,1) & rippletime<=win(iW,2));
    ripplehist{iW} = histcounts(rippletmp,win(iW,1):binsize:win(iW,2));
end

spkhist = cell2mat(spkhist);
ripplehist = cell2mat(ripplehist);

switch iRipple
    case 1
        spkhist(:,ripplehist>0) = [];
    case 2
        spkhist(:,ripplehist==0) = [];
end

nCell = length(spikeTime);
in = false(nCell,nCell);
in(pairIdx(:,1),pairIdx(:,2)) = true;
in = in(:);

x = repmat(unit_id,1,length(unit_id));
y = repmat(unit_id',length(unit_id),1);

pairid = [x(:),y(:)];
pairid = pairid(in,:);

noisecorr = corr(spkhist');
noisecorr = noisecorr(:);
noisecorr = noisecorr(in);

pairid_1 = reshape(pairid(:,1),sum(pairIdx(:,1)),sum(pairIdx(:,2)));
pairid_2 = reshape(pairid(:,2),sum(pairIdx(:,1)),sum(pairIdx(:,2)));
noisecorr = reshape(noisecorr,sum(pairIdx(:,1)),sum(pairIdx(:,2)));
time = length(spkhist)*binsize;
end