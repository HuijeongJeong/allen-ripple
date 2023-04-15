clearvars; clc; close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');

sessionList = session_metric.session_id(session_metric.global_ripple_number>0 &...
    session_metric.immobile_period>100);
nS = length(sessionList);

area = [tag.area.vis6 strcmp(tag.info.structure,'CA1')];
celltype = [tag.celltype.rs tag.celltype.fs];

winSize = 0.05;
winStep = -2:winSize:2;
nWinBin = length(winStep);
resolution = 1;

probeList = {'medial';'lateral'};

% iArea = 2;
% iType = 1;
% parpool(4)
for iS = 1:nS
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_Events.mat'],'pupil_data','running_speed');
    load([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_CA1_ripple','spontaneous_win');
    
%     inVis = ismember(T.unit_id,tag.info.unit_id(area(:,iArea) & celltype(:,iType)));
    timeBin = spontaneous_win(1):winSize:spontaneous_win(2);
    nTimeBin = length(timeBin);
    
    [varData] = cell(1,4);
    
    ripples  = cell(2,1);
    for iProbe = 1:2
        ripples{iProbe} = spontaneous_CA1_ripple.ripple{strcmp(spontaneous_CA1_ripple.relative_location,probeList{iProbe})};
        varData{iProbe} = zeros(nTimeBin-1,nWinBin);
        for iB = 1:nWinBin
            varData{iProbe}(:,iB) = logical(histcounts(ripples{iProbe}(:,1),timeBin-winStep(iB)));
        end
    end
    
    [varData{3},varData{4}] = deal(NaN(nTimeBin-1,1));
    pupil_size = (pupil_data.pupil_height/2.*pupil_data.pupil_width/2*pi)/1000;
    for iB = 1:nTimeBin-1
        iB
        in = running_speed.time>=timeBin(iB) & running_speed.time<timeBin(iB+1);
        varData{3}(iB) = nanmean(running_speed.velocity(in));
        
        in = pupil_data.time>=timeBin(iB) & pupil_data.time<timeBin(iB+1);
        varData{4}(iB) = nanmean(pupil_size(in));
    end
    varData(3:4) = cellfun(@(x) conv(x,fspecial('Gaussian',[1 5*resolution],resolution),'same'),...
        varData(3:4),'UniformOutput',false);
    nVarBin = cellfun(@(x) size(x,2), varData);
    X = cell2mat(varData);
    
    [pval,beta,sem] = deal(cell(size(T,1),4));
    spike = cellfun(@(x) histcounts(x,timeBin),T.spike_time,'UniformOutput',false);
    spikeconv = cellfun(@(x) conv(x,fspecial('Gaussian',[1 5*resolution],resolution),'same'),spike,'UniformOutput',false);
    parfor iC = 1:size(T,1)
        iC
%         [~,~,stats] = glmfit(X,spike{iC}','poisson');
        [~,~,stats] = glmfit(X,spikeconv{iC}','poisson');
        beta(iC,:) = mat2cell(stats.beta(2:end)',1,nVarBin);
        sem(iC,:) = mat2cell(stats.se(2:end)',1,nVarBin);
        pval(iC,:) = mat2cell(stats.p(2:end)',1,nVarBin);
    end
    T1 = cell2table(num2cell(T.unit_id));
    T1.Properties.VariableNames = {'unit_id'};
    
    T2_label = {'b_mrp','b_lrp','b_vel','b_pupil'};
    T2 = cell2table(beta);
    T2.Properties.VariableNames = T2_label;
    
    T3_label = {'se_mrp','se_lrp','se_vel','se_pupil'};
    T3 = cell2table(sem);
    T3.Properties.VariableNames = T3_label;
    
    T4_label = {'p_mrp','p_lrp','p_vel','p_pupil'};
    T4 = cell2table(pval);
    T4.Properties.VariableNames = T4_label;
    
    fr_ripple_glm = [T1,T2, T3, T4];
    save([sdir(sessionList(iS)),'_ripple_modulated.mat'],'fr_ripple_glm','-append');
end

spkRipple = cellfun(@(x) spikeWin(x,[ripples{1}(:,1);ripples{2}(:,1)],[-2 2]),...
    T.spike_time,'UniformOutput',false);
rippleIndex = false(sum(cellfun(@(x) size(x,1),ripples)),2);
rippleIndex(1:size(ripples{1},1),1) = true;
rippleIndex(size(ripples{1},1)+1:end,2) = true;

[xpt,ypt,time,spkhist,spkconv,~,spksem] = cellfun(@(x) rasterPSTH(cellfun(@(y) y*1000,x,'UniformOutput',false),...
    rippleIndex,[-2 2]*1000,winSize*1000,resolution,1),spkRipple,'UniformOutput',false);
time = time{1};

clr = {'r';'b'};
for iC =30:size(T,1)
    axes('Position',axpt(2,2,1,1));
    hold on;
    for iProbe = 1:2
    scatter(xpt{iC}{iProbe}/1000,ypt{iC}{iProbe},2,clr{iProbe},'filled');
    end
    plot([0 0],[0 size(rippleIndex,1)],'k:');
    xlim([-1.5 1.5]);
    ylim([0 size(rippleIndex,1)]);
    title(T.ecephys_structure_acronym(iC));
        
    axes('Position',axpt(2,2,1,2));
    hold on;
    for iProbe = 1:2
       fill([time flip(time)]/1000,spksem{iC}(iProbe,:),clr{iProbe},'EdgeColor','none')
       plot(time/1000,spkconv{iC}(iProbe,:),'Color',clr{iProbe});
    end
    alpha(0.2);
    plot([0 0],get(gca,'YLim'),'k:');
    xlim([-1.5 1.5]);
    
    axes('Position',axpt(2,2,2,1));
    hold on;
    for iProbe = 1:2
    plot(winStep,log10(pval{iC,iProbe}),clr{iProbe});
    if sum(pval{iC,iProbe}(winStep>=-0.25 & winStep<=0.25)<0.01)>2
       text(0,-2.5-(iProbe-1)*0.5,[probeList{iProbe},' modulated']); 
    end
    end
    plot([0 0],[-3 0],'k:');
    plot([-2 2],repmat(log10(0.01),1,2),'k--');
    plot([-2 2],repmat(log10(0.05),1,2),'k--');
    xlim([-1.5 1.5])
    ylim([-5 0]);
    
    
    axes('Position',axpt(2,2,2,2))
    hold on;
    for iProbe = 1:2
        plot(winStep,beta{iC,iProbe},clr{iProbe});
    end
    plot([0 0],[-0.5 0.5],'k:');
    xlim([-1.5 1.5])
    ylim([-0.5 0.5]);
    
    next = input('continue ? (1:yes, 0:no)');
    if next==0
        break
    end
    clf;
end

