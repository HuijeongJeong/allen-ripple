clearvars; clc; close all;

load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
layertable = readtable('D:\OneDrive\1.allen-andermann\layer_info.csv');

areaList = {'VISp';'VISl';'VISrl';'VISal';'VISpm';'VISam'};

[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
session_id = tag.info.session_id(idx);
celltype = [tag.celltype.rs(idx),tag.celltype.fs(idx)];
area = tag.info.structure(idx);
areaidx = cellfun(@(x) find(strcmp(x,areaList)),area,'UniformOutput',false);
areaidx(cellfun(@isempty,areaidx)) = {nan};
areaidx = cell2mat(areaidx);
depth = layertable.cortical_depth(idx);

sessionList = unique(session_id);
nclst =3;
iCT = 1;

ct = cbrewer('qual','Dark2',3);
resolution = 10;

%%
iS = 9;
load([sdir(sessionList(iS)),'_cellTable.mat']);
load([sdir(sessionList(iS)),'_stimulus_evoked.mat'],'fr_gratings');

[in,idx] = ismember(T.unit_id,unit_id.vis(celltype(:,iCT)));
spktime = T.spike_time(in);
clstidx = cluster_idx.vis{nclst-1}(idx(in));
unitarea = areaidx(idx(in));
unitdepth = depth(idx(in));

[~,sortidx] = sortrows([clstidx,unitdepth]);
% [~,sortidx] = sortrows([unitarea,unitdepth]);
clstidx = clstidx(sortidx);
spktime = spktime(sortidx);
unitarea = unitarea(sortidx);
unitdepth = unitdepth(sortidx);
spkhist = fr_gratings.psth.hist(in);
spkhist = spkhist(sortidx);
spkpsth = cell2mat(cellfun(@(x) conv(mean(x),fspecial('Gaussian',...
    [1 5*resolution],resolution),'same'),spkhist,'UniformOutput',false));
%%
ct = cbrewer('qual','Dark2',3);
ct2 = cbrewer('div','RdBu',50);
ct2(ct2<0) = 0;

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5 12]);
h(1) = axes('Position',axpt(8,1,1:7,1));
imagesc(fr_gratings.psth.time,1:sum(in),zscore(spkpsth,[],2));
hold on;
plot([0 0 NaN 2 2],[1 sum(in) NaN 1 sum(in)],'k:');
axis xy
xlim([-0.5 2.5])
% colormap(h(1),ct2);

h(2) = axes('Position',axpt(8,1,8,1));
imagesc(clstidx);
axis xy
colormap(h(2),ct);

set(h,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35);
set(h(2),'YTickLabel',[],'XTick',[]);
set(h(1),'CLim',[-2 2],'XTick',0:2);
ylabel(h(1),'Neuron');
xlabel(h(1),'Time (s)');
cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig3');
print(fHandle,'-depsc','-painters','psth_gratings.ai');


%%

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20 6]); 
h(1) = axes('Position',axpt(1,7,1,1:3)); hold on;
for iunit = 1:length(spktime)
   scatter(spktime{iunit},ones(length(spktime{iunit}),1)*iunit,1,ct(clstidx(iunit),:),'.'); 
end
ylim([1 length(spktime)]);

h(2) = axes('Position',axpt(1,7,1,4:5)); hold on;
spkhist = cell2mat(cellfun(@(x) conv2(histc(x,win(1)-1:bin:win(2)+1)'/bin,...
    fspecial('Gaussian',[1 5*resolution],resolution),'same'),spktime,'UniformOutput',false));
for iclst = 1:3
    m = mean(spkhist(clstidx==iclst,:));
    s = std(spkhist(clstidx==iclst,:))/sqrt(sum(clstidx==iclst));
    fill([win(1)-1:bin:win(2)+1 flip(win(1)-1:bin:win(2)+1)],...
        [m+s flip(m-s)],ct(iclst,:),'EdgeColor','none');
    plot(win(1)-1:bin:win(2)+1,m,'Color',ct(iclst,:));
end

h(3) = axes('Position',axpt(1,7,1,6)); hold on;
plot(running_speed.time,running_speed.velocity_conv,'k');

[~,idx] = histc(pupil_data.time,running_speed.time);
idx(idx==0) = 1;
pupilsize = (pupil_data.pupil_height/2).*(pupil_data.pupil_width/2).*pi;
pupil_running = ~running_speed.immobile(idx);
pupilsize = pupilsize/nanmean(pupilsize(pupil_running &...
    pupil_data.time>=spontaneous_win(1) & pupil_data.time<=spontaneous_win(2)));

h(4) = axes('Position',axpt(1,7,1,7)); hold on;
plot(pupil_data.time,pupilsize,'k');

set(h,'XLim',win,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,...
    'XTick',win(1):10:win(2),'XtickLabel',0:10:diff(win))
set(h(1:3),'XTickLabel',[]);

ylabel(h(1),'Neuron');
ylabel(h(2),'Rate (Hz)');
ylabel(h(3),{'Velocity';'(cm/s)'});
ylabel(h(4),{'Norm.';'pupil'});
xlabel(h(4),'Time (s)');
alpha(h(2),0.2);

cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig2');
print(fHandle,'-depsc','-painters','example_spktrain.ai');
