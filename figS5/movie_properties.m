clc; clearvars; close all;

rng(3);

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
Tfst = readtable('D:\OneDrive\1.allen-andermann\time_to_first_spike.csv');

refList = {'medial';'lateral'};
typeList = {'rs';'fs'};

nclst = 3;

[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
session_id = tag.info.session_id(idx);

unitid = unit_id.vis;
sessionList = unique(session_id);
nS = length(sessionList);

[pval_rf,cluster,area,tfs,rfsize,reliability,selectivity,similarity] = deal(cell(nS,1));
M = 1000;

for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_stimulus_evoked.mat'],'fr_movie','fr_movie_shuffle');
    
    invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis & tag.celltype.rs));
    
    [in,idx] = ismember(T.unit_id,unitid);
    clstidx = NaN(length(T.unit_id),1);
    clstidx(in) = cluster_idx.vis{nclst-1}(idx(in));
    
    [~,idx] = ismember(T.unit_id,tag.info.unit_id);
%     celltype = [tag.celltype.rs(idx), tag.celltype.fs(idx)];
  
%     ctype{iS} = celltype(invis,:);
    cluster{iS} = clstidx(invis);
    area{iS} = T.ecephys_structure_acronym(invis);
    pval_rf{iS} = T.p_value_rf(invis);
        
    [~,idx] = ismember(T.unit_id(invis),Tfst.ecephys_unit_id);
    tfs{iS} = Tfst.time_to_first_spike_fl(idx);
    rfsize{iS} = T.area_rf(invis);
    
    spkhist = cellfun(@(x) cellfun(@(y) histc(y,0:1:30)',x,'UniformOutput',false),...
        fr_movie.spikeTime(invis),'UniformOutput',false);
    spkhist = cellfun(@(x) cell2mat(cellfun(@(y) y(:)',x,'UniformOutput',false)),...
        spkhist,'UniformOutput',false);
    spkhist = cellfun(@(x) x(:,1:end-1),spkhist,'UniformOutput',false);
    spkave = cell2mat(cellfun(@mean,spkhist,'UniformOutput',false));
    
    spkhist_sf = cellfun(@(x) cellfun(@(y) histc(y,0:1:30)',x,'UniformOutput',false),...
        fr_movie_shuffle.spikeTime(invis),'UniformOutput',false);
    spkhist_sf = cellfun(@(x) cell2mat(cellfun(@(y) y(:)',x,'UniformOutput',false)),...
        spkhist_sf,'UniformOutput',false);
    spkhist_sf = cellfun(@(x) x(:,1:end-1),spkhist_sf,'UniformOutput',false);
    
    
    [selectivity{iS}, similarity{iS}, reliability{iS}] = deal(nan(sum(invis),1)); 
    % selectivity R. Quian Quiroga et al, 2007
    % reliability, similarity (vs. shuffled) Rikhye and Sur, 2015
    
    for iC = 1:sum(invis)
       T = min(spkave(iC,:))+[1:M]*(max(spkave(iC,:))-min(spkave(iC,:)))/M;
       R = mean(spkave(iC,:)-T'>0,2);
       A = sum(R)/M;
       selectivity{iS}(iC) = 1-2*A;
       
       reliability{iS}(iC) = nansum(nansum(tril(corrcoef(spkhist{iC}'),-1)))*2/...
           (size(spkhist{iC},1)^2-size(spkhist{iC},1));

       similarity{iS}(iC) = nansum(nansum(corr(spkhist_sf{iC}',spkhist{iC}')))/...
           (size(spkhist_sf{iC},1)*size(spkhist{iC},1));
    end
end
%%

% ctype = cell2mat(ctype);
clusterp = cell2mat(cluster);
clusterp(isnan(clusterp)) = 0;

areap = cat(1,area{:});
% pval_rf = cell2mat(pval_rf);
% tfs = cell2mat(tfs);
% rfsize = cell2mat(rfsize);

% selectivity = cell2mat(selectivity);
% reliability = cell2mat(reliability);
% similarity = cell2mat(similarity);

nCell = 100;
nIter = 100;

%%
close all
sigrf = cellfun(@(x,y) x<0.01&y<2500,pval_rf,rfsize,'UniformOutput',false);

areaLabel = {'V1';'LM';'RL';'AL';'PM';'AM'};
areaList = {'VISp';'VISl';'VISrl';'VISal';'VISpm';'VISam'};
nArea = length(areaList);
datasub = cell(nclst+1,nArea);
ct = [0.6 0.6 0.6;cbrewer('qual','Dark2',3)];

intest = cellfun(@(x,z) x&ismember(z,areaList),sigrf,area,'UniformOutput',false);

group1 = cellfun(@(x,y) x(y),cluster,intest,'UniformOutput',false);
group2 = cellfun(@(x,y) x(y),area,intest,'UniformOutput',false);
group = {cell2mat(group1),categorical(cat(1,group2{:}))};
group{1}(isnan(group{1})) = 0;

for iD = 1:3
    switch iD
        case 1
            data = cell2mat(selectivity);
            binrange = -1:0.1:1;
             ylimit = [0 0.5];
             xticks = -1:1:1;
            datalist = 'Selectivity';
        case 2
            data = cell2mat(reliability);
            binrange = -0.05:0.05:1;
             ylimit = [0 0.5];
             xticks = 0:0.5:1;
            datalist = 'Reliability';
        case 3            
            data = cell2mat(similarity);
            binrange = -0.5:0.05:0.5;
             ylimit = [-0.1 0.1];
             xticks = -0.5:0.5:0.5;
            datalist = 'Similarity';
    end
    
    bincount = cell(nclst+1,1);
    for iClst = 1:nclst+1
        for iIter = 1:nIter
            in = find(clusterp==iClst-1 & cell2mat(sigrf));
            bincount{iClst}(iIter,:) = cumsum(hist(randsample(data(in,:),nCell),binrange));
        end
        for iA = 1:nArea
            in = find(clusterp==iClst-1 & strcmp(areap,areaList{iA}) & cell2mat(sigrf));
            datasub{iClst,iA} = data(in,:);
        end
    end
    
    fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.3 20/6.5]);
    axes('Position',axpt(9,1,1:5,1,[],[0.18 0.05]));
    hold on;
    for iClst = 1:nclst+1
        m = mean(bincount{iClst}/nCell);
        s = sort(bincount{iClst}/nCell);
        s = s([0.05,0.95]*nIter,:);
        fill([binrange flip(binrange)],[s(1,:) flip(s(2,:))],ct(iClst,:),'EdgeColor','none');
        plot(binrange,m,'Color',ct(iClst,:));
    end
    if iD==1
    plot([0 0],[0 1],'k:','LineWidth',0.35);end
    set(gca,'Box','off','TickDir','out','YTick',0:0.2:1,...
        'FontSize',7,'XLim',[min(binrange) max(binrange)],'XTick',xticks)
        
    ylabel('Cumulative fraction');
    xlabel(datalist);
    alpha(0.2);
    
    axes('Position',axpt(9,1,6:9,1,[],[0.18 0.05]));
    hold on;
    s = cellfun(@(x) nanstd(x)/sqrt(length(x)),datasub);
    m = cellfun(@nanmean,datasub);
    for iClst = 1:nclst+1
            errorbar([1:nArea]+0.15*(iClst-1),m(iClst,:),s(iClst,:),'Color',ct(iClst,:),'CapSize',0)
   end
     set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,...
         'XTick',[1:nArea]+0.15,'XTickLabel',areaLabel,'XTickLabelRotation',45,...
         'YTick',ylimit(1):diff(ylimit)/4:ylimit(2),'YLim',ylimit,'XLim',[0.5 6.95]);
    ylabel(datalist);
    
    [panova,~,stat] = anovan(data(cell2mat(intest)),group,'model','interaction','display','off');
    panova
     text(1,diff(ylimit)*0.95+ylimit(1),['P_c_l_u_s_t_e_r = ',num2str(round(panova(1)*1000)/1000)],'FontSize',6);
     text(1,diff(ylimit)*0.87+ylimit(1),['P_a_r_e_a = ',num2str(round(panova(2)*1000)/1000)],'FontSize',6);
    text(1,diff(ylimit)*0.79+ylimit(1),['P_c_X_a = ',num2str(round(panova(3)*1000)/1000)],'FontSize',6);
    
    cd('D:\OneDrive - University of California, San Francisco\figures\allen\figS5');
    print(fHandle,'-depsc','-painters',['movie_',datalist,'.ai']);
end

%% 