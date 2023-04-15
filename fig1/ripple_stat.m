clearvars; clc; close all

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
celltype = [tag.celltype.rs(idx), tag.celltype.fs(idx)];
sessionid = tag.info.session_id(idx);
sessionList = unique(sessionid);
area = tag.info.structure(idx);
nS = length(sessionList);
fraction_global = nan(nS,2);
fraction_sig = nan(nS,2);
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_Events.mat'],'probe');
    load([sdir(sessionList(iS)),'_ripples.mat'],'CA1_ripple_classified','spontaneous_win','spontaneous_anal_win');
    load([sdir(sessionList(iS)),'_cellTable.mat'],'T');
    
    fraction_global(iS,1) = size(CA1_ripple_classified.global_m,1)/...
        (size(CA1_ripple_classified.global_m,1)+size(CA1_ripple_classified.medial,1));
    fraction_global(iS,2) = size(CA1_ripple_classified.global_l,1)/...
        (size(CA1_ripple_classified.global_l,1)+size(CA1_ripple_classified.lateral,1));
    
    fraction_sig(iS,1) = sum(sessionid==sessionList(iS) & celltype(:,1))/...
        sum(tag.info.session_id==sessionList(iS) & tag.celltype.rs & tag.area.vis);
    fraction_sig(iS,2) = sum(sessionid==sessionList(iS) & celltype(:,2))/...
        sum(tag.info.session_id==sessionList(iS) & tag.celltype.fs & tag.area.vis);
    
    fraction_sig(iS,3) = sum(sessionid==sessionList(iS) & celltype(:,1) & ismember(area,{'VISl','VISal'}))/...
        sum(tag.info.session_id==sessionList(iS) & tag.celltype.rs & ismember(tag.info.structure,{'VISl','VISal'}));
    fraction_sig(iS,4) = sum(sessionid==sessionList(iS) & celltype(:,1) & ismember(area,{'VISam','VISpm'}))/...
        sum(tag.info.session_id==sessionList(iS) & tag.celltype.rs & ismember(tag.info.structure,{'VISam','VISpm'}));
end

figure;