clc; clearvars; close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

refList = {'medial';'lateral'};
typeList = {'rs';'fs'};

nclst = 3;

[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
session_id = tag.info.session_id(idx);

unitid = unit_id.vis;
sessionList = unique(session_id);
nS = length(sessionList);

orientationlist = [0,45,90,135];
[clstidx,inanal,prefori] = deal(cell(nS,1));
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_stimulus_evoked.mat'],'beh_gratings','fr_gratings');
    load([sdir(sessionList(iS)),'_Events.mat'],'drifting_gratings','stimulus_condition');
        
     invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis & tag.celltype.rs));
    
    [in,idx] = ismember(T.unit_id(invis),unitid);
    clstidx{iS} = NaN(sum(invis),1);
    clstidx{iS}(in) = cluster_idx.vis{nclst-1}(idx(in));
    
    trialinfo = drifting_gratings;
    frinfo = fr_gratings;
    behinfo = beh_gratings;
    
%     [~,idx] = ismember(trialinfo.stimulus_condition_id,stimulus_condition.stimulus_condition_id);
%     orientation = cellfun(@str2num,stimulus_condition.orientation(idx));

    %%
    inbasetime = fr_gratings.psth.time>=-0.5 & fr_gratings.psth.time<0;
    instimtime = fr_gratings.psth.time>=0 & fr_gratings.psth.time<=2;
    %%
    frbase = cellfun(@(x,y) mean(x(drifting_gratings.stimulus_condition_id==y,...
        inbasetime),2),fr_gratings.psth.hist(invis),...
        num2cell(fr_gratings.pref_condition_id(invis)),'UniformOutput',false);
    frpref = cellfun(@(x,y) mean(x(drifting_gratings.stimulus_condition_id==y,...
        instimtime),2),fr_gratings.psth.hist(invis),...
        num2cell(fr_gratings.pref_condition_id(invis)),'UniformOutput',false);
    prefori{iS} = cellfun(@(x) str2double(stimulus_condition.orientation(stimulus_condition.stimulus_condition_id==x)),...
        num2cell(fr_gratings.pref_condition_id(invis)));
    
    [~,p_responsiveness] = cellfun(@(x,y) ttest(x,y),frbase,frpref);
    sig_dg = p_responsiveness<0.05;
    sig_rf = T.p_value_rf(invis)<0.01 & T.area_rf(invis)<2500;
    inanal{iS} = sig_dg & sig_rf;
    %%
%     avefr = nan(sum(invis),4);
%     for io = 1:4
%         avefr(:,io) = cellfun(@(x) mean(mean(x(orientation==orientationlist(io),instimtime),2)),...
%             fr_gratings.psth.hist(invis));
%     end
%     [~,prefori{iS}] = max(avefr,[],2);
end

%%

clstidx_t = cell2mat(clstidx);
prefori_t = cell2mat(prefori);
inanal_t = cell2mat(inanal);
data = nan(4,3);
for iclst = 1:3
   for io = 1:4
      data(io,iclst) = sum(clstidx_t==iclst&prefori_t==orientationlist(io)&inanal_t); 
   end
end
prob = data./repmat(sum(data,1),4,1);
%%
c_cluster = cbrewer('qual','Dark2',7);
c_cluster = c_cluster(4:7,:);
close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 4]);
hold on;
for iclst = 1:3
    b = bar(iclst,prob(:,iclst),'stacked', 'FaceColor','flat');
    for io = 1:4
        b(io).CData = c_cluster(io,:);
        text(iclst-0.3,mean([sum(prob(1:io-1,iclst)),sum(prob(1:io,iclst))]),...
            [num2str(data(io,iclst)),newline,'(',...
            num2str(round(prob(io,iclst)*1000)/10),'%)'],'FontSize',5);
    end
end
set(gca,'XLim',[0.2 3.8],'Box','off','TickDir','out','XTick',1:3,...
    'XTickLabel',{'iAct','dAct','Inh'},'XTickLabelRotation',45,'YTick',0:0.5:1,'YTickLabel',[0,50,100]);
ylabel('% neuron')
%print(fHandle,'-dtiff','-r600','D:\OneDrive - UCSF\data\allen\revision\clustering_osi.tif');
print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\figures\2.allen\revision\FigS_visualresponse\clustering_osi.ai');
