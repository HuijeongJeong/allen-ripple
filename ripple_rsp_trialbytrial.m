clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

winbase = [-3, -2];
win = [-3, 3];
binsize = 0.2;
rippletype = {'medial','lateral','global'};
clstname = {'Thal','iAct','dAct','Inh'};
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_cellTable.mat'],'T');
    load([sdir(sessionList(iS)),'_ripples.mat'],'CA1_ripple_classified');
    
    targets = tag.info.unit_id((ismember(tag.info.unit_id,unit_id.vis)&tag.celltype.rs) |...
        tag.area.thalamus);
    in = ismember(T.unit_id,targets);
    
    unitid = T.unit_id(in);
    spiketime = T.spike_time(in);
    
    cidx = zeros(sum(in),1);
    [in,idx] = ismember(unitid,unit_id.vis);
    cidx(in) = cluster_idx.vis{2}(idx(in));
    
    %%
    fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 10 10]);
    for iRp = 1:3
        ripple = CA1_ripple_classified.(rippletype{iRp})(:,1);
        spkbaseline = cellfun(@(y) cellfun(@(x) sum(y>=x+winbase(1) & y<=x+winbase(2)),...
            num2cell(ripple)),spiketime,'UniformOutput',false);
        spktest = cellfun(@(y) cell2mat(cellfun(@(x) histc(y,x+[win(1):binsize:win(2)])',...
            num2cell(ripple),'UniformOutput',false)),spiketime,'UniformOutput',false);
        normspk = cellfun(@(x,y) y(:,1:end-1)-repmat(x/(diff(winbase)/binsize),1,diff(win)/binsize),...
            spkbaseline,spktest,'UniformOutput',false);
        normspk = cat(3,normspk{:});
        %%
        for iclst = 1:4
            axes('Position',axpt(4,3,iclst,iRp))
            avedeltaspk = mean(normspk(:,:,cidx==iclst-1),3);
            imagesc(win(1):binsize:win(2),1:length(ripple),avedeltaspk);
            hold on;
            plot([0 0],[0 length(ripple)],'k:');
            set(gca,'CLim',[-1.5 1.5],'FontSize',8);
            if iRp==1
                title(clstname{iclst});
            end
            if iclst==1
                ylabel(rippletype{iRp});
            end
        end
    end
    print(fHandle,'-dtiff','-r600',[sdir(sessionList(iS)),'_ripplersp.tif'])
end
