clearvars; clc; close all;

load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
layertable = readtable('D:\OneDrive\1.allen-andermann\layer_info.csv');

%%
nclst =3;

unitid = tag.info.unit_id(tag.area.vis & tag.celltype.rs);
sessionid = tag.info.session_id(tag.area.vis & tag.celltype.rs);
% celltype = [tag.celltype.rs(tag.area.vis),tag.celltype.fs(tag.area.vis)];
% depth = layertable.cortical_depth(tag.area.vis);
[in,idx] = ismember(unitid,unit_id.vis);
clusteridx = zeros(length(unitid),1);
clusteridx(in) = cluster_idx.vis{nclst-1}(idx(in));

sessionList = unique(sessionid);
ct = [0.6 0.6 0.6; cbrewer('qual','Dark2',3)];

%%
iS = 9;
load([sdir(sessionList(iS)),'_cellTable.mat']);
load([sdir(sessionList(iS)),'_Events'],'running_speed','pupil_data');
load([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_win','spontaneous_anal_win',...
    'CA1_ripple_classified','filtered_lfp');

%%
[in,idx] = ismember(T.unit_id,unitid);
spktime = T.spike_time(in);
firingrate_spon = cellfun(@(x) sum(x>=spontaneous_win(1) & x<=spontaneous_win(2))/...
    diff(spontaneous_win),T.spike_time(in));

clstidx = clusteridx(idx(in));
[~,sortidx] = sortrows([clstidx,firingrate_spon]);
clstidx = clstidx(sortidx);
spktime = spktime(sortidx);

out = randsample(find(clstidx==0),sum(clstidx==0)-30);
clstidx(out) = [];
spktime(out) = [];

for iC = 1:3
    if sum(clstidx==iC)>30
        out = find(find(clstidx==iC),sum(clstidx==iC)-30,'first');
     end
   clstidx(out) = [];
    spktime(out) = [];
end

%%
win = [4417 4457]+3;
% win = [5815 5865];
bin = 0.01;
resolution = 10;
rippletype = {'medial','lateral','global'};
cripple = {'r';'b';'k'};

close all

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20 6]); 
h(1) = axes('Position',axpt(1,7,1,1:4)); hold on;
for iunit = 1:length(spktime)
%     if clstidx(iunit)==0
%         continue;
%     end
   scatter(spktime{iunit},ones(length(spktime{iunit}),1)*iunit,...
       1,ct(clstidx(iunit)+1,:),'.'); 
end


ylim([1 length(spktime)+1]);

%%
h(2) = axes('Position',axpt(1,7,1,5)); hold on;
% spkhist = cell2mat(cellfun(@(x) conv2(zscore(histc(x,win(1)-1:bin:win(2)+1)'/bin),...
%     fspecial('Gaussian',[1 5*resolution],resolution),'same'),spktime,'UniformOutput',false));
% for iclst = 1:4
%     m = mean(spkhist(clstidx==iclst-1,:));
%     s = std(spkhist(clstidx==iclst-1,:))/sqrt(sum(clstidx==iclst-1));
%     fill([win(1)-1:bin:win(2)+1 flip(win(1)-1:bin:win(2)+1)],...
%         [m+s flip(m-s)],ct(iclst,:),'EdgeColor','none');
%     plot(win(1)-1:bin:win(2)+1,m,'Color',ct(iclst,:));
% end
% for iRp = 1:3
%    rippleonset = CA1_ripple_classified.(rippletype{iRp})(:,1);
%    inripple = rippleonset>=win(1) & rippleonset<=win(2);
%    x = [repmat(rippleonset(inripple),1,2),nan(sum(inripple),1)]';
%    y = repmat([-0.3;0.6; NaN],1,sum(inripple));
%    plot(x(:),y(:),cripple{iRp},'LineStyle',':');
% end
% ylim([-0.3 0.6]);
for iRp = 1:3
    if iRp<3
nss = mean(NormalizedSquaredSignal_HJ(...
    [filtered_lfp.time{iRp}',filtered_lfp.lfp{iRp}]),2);
time = filtered_lfp.time{iRp};
plot(time,nss,cripple{iRp});
    end
rippleonset = CA1_ripple_classified.(rippletype{iRp})(:,1);
scatter(rippleonset,ones(length(rippleonset),1)*29,2,cripple{iRp},'filled');
end

%%
h(3) = axes('Position',axpt(1,7,1,6)); hold on;

plot(running_speed.time,running_speed.velocity,'k');


h(4) = axes('Position',axpt(1,7,1,7)); hold on;

[~,idx] = histc(pupil_data.time,running_speed.time);
idx(idx==0) = 1;
pupilsize = (pupil_data.pupil_height/2).*(pupil_data.pupil_width/2).*pi;
pupil_running = ~running_speed.immobile(idx);
pupilsize = pupilsize/nanmean(pupilsize(pupil_running &...
    pupil_data.time>=spontaneous_win(1) & pupil_data.time<=spontaneous_win(2)));

% yyaxis right
plot(pupil_data.time,pupilsize,'k');
set(h,'XLim',win);
set(h,'XLim',win,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,...
    'XTick',win(1):10:win(2),'XtickLabel',0:10:diff(win))
set(h(1:3),'XTickLabel',[]);

ylabel(h(1),'Neuron');
% ylabel(h(2),'Norm. firing rate (z-score)');
ylabel(h(3),{'Norm. ripple';'band power'});
ylabel(h(4),{'Norm.';'pupil'});
% yyaxis left
ylabel(h(3),{'Velocity';'(cm/s)'});
xlabel(h(4),'Time (s)');
alpha(h(2),0.2);
cd('D:\OneDrive - UCSF\figures\allen\fig2');
% cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig2');
print(fHandle,'-depsc','-painters','example_spktrain.ai');
