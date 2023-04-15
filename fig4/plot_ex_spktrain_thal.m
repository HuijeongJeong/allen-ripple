clearvars; clc; close all;

load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
layertable = readtable('D:\OneDrive\1.allen-andermann\layer_info.csv');

[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
session_id = tag.info.session_id(idx);
celltype = [tag.celltype.rs(idx),tag.celltype.fs(idx)];
depth = layertable.cortical_depth(idx);

sessionList = unique(session_id);
nclst =3;
iCT = 1;

ct = cbrewer('qual','Dark2',3);

%%
iS = 9;
load([sdir(sessionList(iS)),'_cellTable.mat']);
load([sdir(sessionList(iS)),'_Events'],'running_speed','pupil_data');
load([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_win','spontaneous_anal_win');

[in,idx] = ismember(T.unit_id,unit_id.vis(celltype(:,iCT)));
spktime = T.spike_time(in);
clstidx = cluster_idx.vis{nclst-1}(idx(in));

inthal = ismember(T.unit_id,tag.info.unit_id(tag.area.thalamus));
spktime_thal = T.spike_time(inthal);
area_thal = T.ecephys_structure_acronym(inthal);
[area_thal,sortidx] = sort(area_thal);
spktime_thal = spktime_thal(sortidx);

%%
win = [4440 4460];
% win = [4445 4465];
% win = [5815 5865];
bin = 0.01;
resolution = 10;

close all

ct = [1 0 0; 1 0.6 0.6; 0 0 0];

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 10 6]); 
h(1) = axes('Position',axpt(1,7,1,1:3)); hold on;
for iunit = 1:sum(inthal)
    if strcmp(area_thal{iunit},'LGd')
       i = 1;
    elseif strcmp(area_thal{iunit},'LP')
       i = 2;
    else
       i = 3;
    end
   scatter(spktime_thal{iunit},ones(length(spktime_thal{iunit}),1)*iunit,1,ct(i,:),'.'); 
end
ylim([1 sum(inthal)]);
%%
h(2) = axes('Position',axpt(1,7,1,4:5)); hold on;
spkhist = cell2mat(cellfun(@(x) conv2(histc(x,win(1)-1:bin:win(2)+1)'/bin,...
    fspecial('Gaussian',[1 5*resolution],resolution),'same'),spktime_thal,'UniformOutput',false));
for ithal = 1:3
    if ithal==1
       in = cellfun(@(x) strcmp(x,'LGd'),area_thal); 
    elseif ithal==2
       in = cellfun(@(x) strcmp(x,'LP'),area_thal);  
    else
       in = cellfun(@(x) ~contains(x,{'LGd','LP'}),area_thal); 
    end
    
    m = mean(spkhist(in,:));
    s = std(spkhist(in,:))/sqrt(sum(in));
    fill([win(1)-1:bin:win(2)+1 flip(win(1)-1:bin:win(2)+1)],...
        [m+s flip(m-s)],ct(ithal,:),'EdgeColor','none');
    plot(win(1)-1:bin:win(2)+1,m,'Color',ct(ithal,:));
end

%%
ct = cbrewer('qual','Dark2',3);
h(3) = axes('Position',axpt(1,7,1,6:7)); hold on;
spkhist = cell2mat(cellfun(@(x) conv2(histc(x,win(1)-1:bin:win(2)+1)'/bin,...
    fspecial('Gaussian',[1 5*resolution],resolution),'same'),spktime,'UniformOutput',false));
for iclst = 1:3
    m = mean(spkhist(clstidx==iclst,:));
    s = std(spkhist(clstidx==iclst,:))/sqrt(sum(clstidx==iclst));
    fill([win(1)-1:bin:win(2)+1 flip(win(1)-1:bin:win(2)+1)],...
        [m+s flip(m-s)],ct(iclst,:),'EdgeColor','none');
    plot(win(1)-1:bin:win(2)+1,m,'Color',ct(iclst,:));
end

%%

set(h,'XLim',win,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,...
    'XTick',win(1):10:win(2),'XtickLabel',0:10:diff(win))
set(h(1:2),'XTickLabel',[]);

ylabel(h(1),'Neuron');
ylabel(h(2),'Rate (Hz)');
ylabel(h(3),'Rate (Hz)');
xlabel(h(3),'Time (s)');
alpha(h(2),0.2);
alpha(h(3),0.2);

cd('D:\OneDrive - University of California, San Francisco\figures\2.allen\fig4');
print(fHandle,'-depsc','-painters','example_spktrain_thal2.ai');
