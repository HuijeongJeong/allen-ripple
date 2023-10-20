clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

win = [-5 5];
binsize = 0.01;
resolution = 10;
mininterval = 0.5;

[cidx,celltype] = deal(cell(nS,1));
iri = cell(nS,5);
% ripplehist = cell(3,4);
data = cell(nS,2,2);
rippletype = {'medial','lateral','global','global_m','global_l'};
numripple = nan(nS,4);

for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_ripple_modulated.mat'],'fr_ripple');
    load([sdir(sessionList(iS)),'_cellTable.mat'],'T');
    load([sdir(sessionList(iS)),'_ripples.mat'],'CA1_ripple_classified');
    
    [in,idx] = ismember(T.unit_id,unit_id.vis);
    cidx{iS} = cluster_idx.vis{2}(idx(in));
    [~,idx] = ismember(T.unit_id(in),tag.info.unit_id);
    celltype{iS} = [tag.celltype.rs(idx), tag.celltype.fs(idx)];
    
    rippletime = [CA1_ripple_classified.medial(:,1);...
        CA1_ripple_classified.lateral(:,1);
        CA1_ripple_classified.global(:,1)];
    rippletypeidx = [ones(size(CA1_ripple_classified.medial,1),1);...
        ones(size(CA1_ripple_classified.lateral,1),1)*2;
        ones(size(CA1_ripple_classified.global,1),1)*3];
    ripple = [rippletime,rippletypeidx];
    
    [ripple,sortidx] = sortrows(ripple,1);
    
    %%
    ripple(1:round(size(ripple,1)/4),3) = 1; % first quadrant
    ripple(end-round(size(ripple,1)/4)+1:end,3) = 2; % last quadrant
   
    %%
    % number of total riples & isolated ripples    
    for iw = 1:2
        spkhist = cellfun(@(y) cell2mat(cellfun(@(x) histc(y,x+[win(1):binsize:win(2)])'*(1/binsize),...
            num2cell(ripple(ripple(:,3)==iw,1)),'UniformOutput',false)),...
            T.spike_time(in),'UniformOutput',false);
        spkhist = cellfun(@(x) x(:,1:end-1),spkhist,'UniformOutput',false);
        
        for irp = 1:2
            spkave = conv2(cell2mat(cellfun(@(x) nanmean(x(ripple(ripple(:,3)==iw,2)==irp,:)),...
                spkhist,'UniformOutput',false)),fspecial('Gaussian',[1 5*resolution],resolution),'same');
            data{iS,iw,irp} = spkave;
        end
        
    end
end

%%
close all
clr = cbrewer('qual','Dark2',3);

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 15 8]);
intime = fr_ripple.time>=-1.5 & fr_ripple.time<=1.5;
time = fr_ripple.time(intime);
for i = 1:2
    clstidx = cell2mat(cellfun(@(x,y) x(y(:,1)),cidx,celltype,'UniformOutput',false));
    dataz = zscore(cell2mat(cellfun(@(x,y) x(y(:,1),intime),...
        data(:,:,i),repmat(celltype,1,2),'UniformOutput',false)),[],2);
    [~,score] = pca(dataz);
    if i==1
        [~,sortidx] = sortrows([clstidx,...
            max(abs(dataz(:,[100:200,400:500])),[],2)],{'descend','ascend'});
    end
    for iRp = 1:2
        %         axes('Position',axpt(7,1,[1:3]+(iRp-1)*3,1,axpt(2,3,i,1:2)));
        axes('Position',axpt(2,1,iRp,1,axpt(15,1,[1:7]+(i-1)*7,1)));
        hold on;
        imagesc(time,1:length(sortidx),dataz(sortidx,[1:300]+(iRp-1)*300));
        plot([0 0],[0.5 length(sortidx)+0.5],'k:');
        axis xy
        axis tight
        set(gca,'CLim',[-2 2],'XTick',-1:1,'YTick',0:400:1600,'Box','off',...
            'TickDir','out','FontSize',8)
        xlabel('Time (s)');
        if iRp==2
            set(gca,'YTickLabel',[]);
        else
            if i==1
                ylabel('Cortical neuron #');
            else
                set(gca,'YTickLabel',[]);
            end
        end
    end
    
    
end
h = axes('Position',axpt(15,1,15,1));
imagesc(clstidx(sortidx));
axis xy
colormap(h,clr);
set(gca,'Box','off','TickDir','out','XTick',[],'YTick',0:400:1600,'YTickLabel',[]);
print(fHandle,'-dtiff','-r600','D:\OneDrive - UCSF\data\allen\revision\psth_firstlast.tif')
