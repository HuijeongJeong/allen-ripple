clc; clearvars; close all;

rng(3);

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
Tfst = readtable('D:\OneDrive\1.allen-andermann\time_to_first_spike.csv');
layertable = readtable('D:\OneDrive\1.allen-andermann\layer_info.csv');

refList = {'medial';'lateral'};
typeList = {'rs';'fs'};

nclst = 3;

[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
session_id = tag.info.session_id(idx);

unitid = unit_id.vis;
sessionList = unique(session_id);
nS = length(sessionList);

resolution = 5;

%%
% [pval_rf,ctype,cluster,area,tfs,rfsize,modidx,fr_dg] = deal(cell(nS,1));

[cluster,area,fraveconvz,responsiveness,selectivity,...
    p_responsive,modindex,pval_rf,area_rf,latency,layer,depth] = deal(cell(nS,1));
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_stimulus_evoked.mat'],'fr_gratings');
    load([sdir(sessionList(iS)),'_Events.mat'],'drifting_gratings','stimulus_condition');
    
    %%
    invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis));
    [in,idx] = ismember(T.unit_id,unitid);
    clstidx = zeros(length(T.unit_id),1);
    clstidx(in) = cluster_idx.vis{nclst-1}(idx(in));
    
    [~,idx] = ismember(T.unit_id,tag.info.unit_id);
    celltype = [tag.celltype.rs(idx), tag.celltype.fs(idx)];
    
    in = invis & celltype(:,1);
    
    cluster{iS} = clstidx(in);
    area{iS} = T.ecephys_structure_acronym(in);
    modindex{iS} = T.mod_idx_dg(in);
    pval_rf{iS} = T.p_value_rf(in);
    area_rf{iS} = T.area_rf(in);
    
    %% 
    [~,idx] = ismember(T.unit_id(in),layertable.ecephys_unit_id);
    layer{iS} = layertable.cortical_layer(idx);
    depth{iS} = layertable.cortical_depth(idx);

    %     fr_dg{iS} = T.firing_rate_dg(in);
    
    %%
    movtime = movmean(fr_gratings.psth.time,5,'EndPoints','discard');
    frhist = cellfun(@(x,y) movmean(x(drifting_gratings.stimulus_condition_id==y,:),5,2,'Endpoints','discard'),...
        fr_gratings.psth.hist(in),num2cell(fr_gratings.pref_condition_id(in)),'UniformOutput',false);
    
    
    cellfun(@(x) ttest(repmat(x(:,9),1,size(x,2)),x),frhist,'UniformOutput',false);
    
    frave = cellfun(@(x,y) mean(x(drifting_gratings.stimulus_condition_id==y,:)),...
        fr_gratings.psth.hist(in),num2cell(fr_gratings.pref_condition_id(in)),'UniformOutput',false);
    frbase_pref = cellfun(@(x,y) mean(x(drifting_gratings.stimulus_condition_id==y,...
        fr_gratings.psth.time>=-0.5 & fr_gratings.psth.time<0),2),...
        fr_gratings.psth.hist(in),num2cell(fr_gratings.pref_condition_id(in)),'UniformOutput',false);
    fraveconv = cellfun(@(x) conv(x,fspecial('Gaussian',[1 5*resolution],resolution),'same'),...
        frave,'UniformOutput',false);
    fraveconvz{iS} = cell2mat(cellfun(@(x,y) (x-mean(y))/std(y),fraveconv,frbase_pref,'UniformOutput',false));
    
    %%
    [~,idx] = ismember(fr_gratings.pref_condition_id(in),stimulus_condition.stimulus_condition_id);
    pref_ori = cellfun(@str2num,stimulus_condition.orientation(idx));
    orth_ori = pref_ori-90;
    orth_ori(orth_ori<0) = 180+orth_ori(orth_ori<0);
    pref_contrast = cellfun(@str2num,stimulus_condition.contrast(idx));
    
    indg = stimulus_condition.stimulus_name=="drifting_gratings_75_repeats";
    condition_id = stimulus_condition.stimulus_condition_id(indg);
    orth_condition_id = condition_id(cellfun(@(x,y) find(cellfun(@str2num,stimulus_condition.orientation(indg))==x &...
        cellfun(@str2num,stimulus_condition.contrast(indg))==y), num2cell(orth_ori),num2cell(pref_contrast)));
    %%
    
    frpref = cellfun(@(x,y) mean(x(drifting_gratings.stimulus_condition_id==y,...
        fr_gratings.psth.time>=0 & fr_gratings.psth.time<=2),2),...
        fr_gratings.psth.hist(in),num2cell(fr_gratings.pref_condition_id(in)),'UniformOutput',false);
    frorth = cellfun(@(x,y) mean(x(drifting_gratings.stimulus_condition_id==y,...
        fr_gratings.psth.time>=0 & fr_gratings.psth.time<=2),2),...
        fr_gratings.psth.hist(in),num2cell(orth_condition_id),'UniformOutput',false);
    frbase_orth = cellfun(@(x,y) mean(x(drifting_gratings.stimulus_condition_id==y,...
        fr_gratings.psth.time>=-0.5 & fr_gratings.psth.time<0),2),...
        fr_gratings.psth.hist(in),num2cell(orth_condition_id),'UniformOutput',false);
    
    responsiveness{iS} = (cellfun(@mean,frpref)-cellfun(@mean,frbase_pref))./...
        sqrt(cellfun(@std,frbase_pref).^2+cellfun(@std,frpref).^2);
    selectivity{iS} = (cellfun(@(x,y) mean(x-y),frpref,frbase_pref)-...
        cellfun(@(x,y) mean(x-y),frorth,frbase_orth))./...
        sqrt(cellfun(@(x,y) std(x-y),frpref,frbase_pref).^2+...
        cellfun(@(x,y) std(x-y),frorth,frbase_orth).^2);
    [~,p_responsive{iS}] = cellfun(@(x,y) ttest(x,y),frpref,frbase_pref);
    
    %%
    maxrsp = max(fraveconvz{iS}(:,fr_gratings.psth.time>=0 & fr_gratings.psth.time<=2),[],2);
    minrsp = min(fraveconvz{iS}(:,fr_gratings.psth.time>=0 & fr_gratings.psth.time<=2),[],2);
    frconvtemp = fraveconvz{iS}(:,fr_gratings.psth.time>=0 & fr_gratings.psth.time<=2);
    timetemp = fr_gratings.psth.time(fr_gratings.psth.time>=0 & fr_gratings.psth.time<=2);
    [~,maxidx] = cellfun(@(x,y) find(x>=y*0.5,1,'first'),...
        mat2cell(frconvtemp,ones(size(frconvtemp,1),1),size(frconvtemp,2)),...
        num2cell(maxrsp),'UniformOutput',false);
    [~,minidx] = cellfun(@(x,y) find(x<=y*0.5,1,'first'),...
        mat2cell(frconvtemp,ones(size(frconvtemp,1),1),size(frconvtemp,2)),...
        num2cell(minrsp),'UniformOutput',false);
    latency{iS} = NaN(size(frconvtemp,1),1);
    latency{iS}(responsiveness{iS}<0) = timetemp(cell2mat(minidx(responsiveness{iS}<0)));
    latency{iS}(responsiveness{iS}>0) = timetemp(cell2mat(maxidx(responsiveness{iS}>0)));
end

%%
ct = [0.6 0.6 0.6; cbrewer('qual','Dark2',3)];
sigrf = cellfun(@(x,y) x<0.01 & y<2500,pval_rf,area_rf,'UniformOutput',false);

%% significant receptive field fraction & receptive field size
layerdist = cell2mat(cellfun(@(x,y) histc(x(y>0),1:6)',layer,cluster,'UniformOutput',false));
sum(layerdist)/sum(layerdist(:))
%%
close all

nIter = 100;
nCell = 100;
layerfr = [0.15,0.15,0.45,0.25];
layerlist = [2,4,5,6];

clustert = cell2mat(cluster);
sigrft = cell2mat(sigrf);
layert = cell2mat(layer);

sigrf_clst = nan(nIter,4);
sigrf_clst_layer = nan(nIter,4);
for iIter = 1:nIter
   for iClst = 1:4
      sigrf_clst(iIter,iClst) = mean(sigrft(randsample(find(clustert==iClst-1),nCell)));
      sigrf_clst_layer_temp = [];
      for iLayer = 1:4
          sigrf_clst_layer_temp = [sigrf_clst_layer_temp;...
              sigrft(randsample(find(clustert==iClst-1 &...
          layert==layerlist(iLayer)),nCell*layerfr(iLayer)))];
      end
      sigrf_clst_layer(iIter,iClst) = mean(sigrf_clst_layer_temp);
   end
end

sigrf_clst = sort(sigrf_clst,1);
sigrf_clst_layer = sort(sigrf_clst_layer,1);
dir = 'D:\OneDrive - University of California, San Francisco\figures\allen\fig3';
for i = 1:2
   fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 3.8]);
   hold on;
   x = [repmat(1:4,2,1);nan(1,4)];
   if i==1
      y = [sigrf_clst(0.05*nIter,:);sigrf_clst(0.95*nIter,:);nan(1,4)];
      m = mean(sigrf_clst);
      
   else
      y = [sigrf_clst_layer(0.05*nIter,:);sigrf_clst_layer(0.95*nIter,:);nan(1,4)];
      m = mean(sigrf_clst_layer);
   end
   plot(x(:),y(:),'k');
   scatter(1:4,m,2,'k','filled');
   xlim([0.5 4.5])
   ylim([0 1])
   ylabel('% siginificant receptive field');
   set(gca,'XTick',1:4,'XTickLabel',{'Nomod';'iAct';'dAct';'Inh'},'XTickLabelRotation',45,...
       'Box','off','TickDir','out','FontSize',8,'YTick',0:0.2:1,'YTickLabel',0:20:100);
   if i==1
   print(fHandle,'-depsc','-painters',[dir,'\significant_receptive_fields.ai']);
   else
   print(fHandle,'-depsc','-painters',[dir,'\significant_receptive_fields_layermatched.ai']);
   end
end


%%


binrange = 0:100:2500;
bincount = NaN(4,length(binrange));
nsigrf = NaN(4,1);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.275 20/6.5]);
h(1) = axes('Position',axpt(5,1,1:2,1,[],[0.2 0.05])); hold on;
h(2) = axes('Position',axpt(5,1,3:5,1,[],[0.2 0.05])); hold on;
for iClst = 1:4
    sigrf_clst = cellfun(@(x,y) x(y==iClst-1),sigrf,cluster,'UniformOutput',false);
    arearf_clst = cell2mat(cellfun(@(x,y) x(y),area_rf,sigrf_clst,'UniformOutput',false));
    
    nsigrf(iClst) = mean(cell2mat(sigrf_clst));
    
    scatter(iClst,nsigrf(iClst),4,ct(iClst,:),'filled','Parent',h(1));
%     bincount(iClst,:) = histc(cell2mat(arearf_clst),binrange);
    for iIter = 1:nIter
       bincount(iIter,:) = histc(randsample(arearf_clst,nCell),binrange);
    end
    m = mean(cumsum(bincount,2)./repmat(sum(bincount,2),1,length(binrange)));
    s = sort(cumsum(bincount,2)./repmat(sum(bincount,2),1,length(binrange)));
    s = s([0.05 0.95]*nIter,:);
    fill([binrange,flip(binrange)],[s(1,:) flip(s(2,:))],ct(iClst,:),'EdgeColor','none','Parent',h(2));
    plot(binrange,m,'Color',ct(iClst,:),'Parent',h(2));
%     plot(binrange,cumsum(bincount(iClst,:))/sum(bincount(iClst,:)),'Color',ct(iClst,:),'Parent',h(2));
end
set(h,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'YTick',0:0.2:1);
set(h(1),'XLim',[0.5 4.5],'YLim',[0 0.8],'XTick',1:4,'XTickLabel',{'Nomod','iAct','dAct','Inh'},...
    'XTickLabelRotation',45,'YTickLabel',0:20:100);
set(h(2),'XLim',[0 2500]);
ylabel(h(1),'% significant receptive field');
ylabel(h(2),'Cumulative fraction');
xlabel(h(2),'Receptive field area (deg^2)');
set(h(2),'XTick',0:1000:2000);
cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig3');
% print(fHandle,'-depsc','-painters','receptive_field.ai');

%% fraction responsive to gratings
close all

responsivetype = cellfun(@(x,y,z) [x>=0.05, x<0.05&y>0, x<0.05&y<0] & z,...
    p_responsive,responsiveness,sigrf,'UniformOutput',false);

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.1 20/6.5]);
h(1) = axes('Position',axpt(5,1,2:5,1,[],[0.2 0.05])); hold on;
hold on;
for iClst = 1:4
    fon(iClst,:) = sum(cell2mat(cellfun(@(x,y) x(y==iClst-1,:),responsivetype,cluster,...
        'UniformOutput',false)));
end
% plot(1:4,1-fon(:,1)./sum(fon,2),'k');
plot(1:4,fon(:,2)./sum(fon,2),'r');
plot(1:4,fon(:,3)./sum(fon,2),'b');

ylabel('% responsive to gratings');
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',1:4,...
    'XTickLabel',{'Nomod','iAct','dAct','Inh'},'XTickLabelRotation',45,...
    'XLim',[0.5 4.5],'YLim',[0 0.8],'YTick',0:0.2:1,'YTickLabel',0:20:100);
print(fHandle,'-depsc','-painters','responsive_to_dg.ai');


%% heatmap for gratings
ct = [0.6 0.6 0.6; cbrewer('qual','Dark2',3)];
ct2 = [1 1 1; 1 0 0; 0 0 1];

% close all
% sortdata = cell2mat(cellfun(@(x,y,z,w) [x(strcmp(z,'VISp')),y(strcmp(z,'VISp'))],...
%     cluster,responsiveness,area,p_responsive,'UniformOutput',false));
% [sortdata,sortidx] = sortrows(sortdata);
% data = cell2mat(cellfun(@(x,y,w) x(strcmp(y,'VISp') & w<0.01,:),fraveconvz,area,p_responsive,'UniformOutput',false));


close all
sortdata = cell2mat(cellfun(@(x,y,z,w) [x(z),y(z),double(w(z)<0.05)],cluster,responsiveness,...
    sigrf,p_responsive,'UniformOutput',false));
sortdata(sortdata(:,2)<0 & sortdata(:,3),3) = 2;
[sortdata,sortidx] = sortrows(sortdata);
data = cell2mat(cellfun(@(x,y) x(y,:),fraveconvz,sigrf,'UniformOutput',false));

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 10 7]);
h(1) = axes('Position',axpt(10,1,1:8,1,axpt(2,1,1,1,[],[0.1 0.05])));
imagesc(fr_gratings.psth.time,1:sum(sortdata(:,1)==0),data(sortidx(sortdata(:,1)==0),:));
hold on;
plot([0 0 NaN 2 2 ],[1 sum(sortdata(:,1)==0) NaN 1 sum(sortdata(:,1)==0)],'k:');
xlabel('Time from grating onset (s)');
ylabel('Cortical neuron #');
h(2) = axes('Position',axpt(10,1,9,1,axpt(2,1,1,1,[],[0.1 0.05])));
imagesc(ones(sum(sortdata(:,1)==0),1));
colormap(h(2),[0.6 0.6 0.6]);
h(5) = axes('Position',axpt(10,1,10,1,axpt(2,1,1,1,[],[0.1 0.05])));
imagesc(sortdata(sortdata(:,1)==0,3));
colormap(h(5),ct2);

h(3) = axes('Position',axpt(10,1,1:8,1,axpt(2,1,2,1,[],[0.1 0.05])));
imagesc(fr_gratings.psth.time,1:sum(sortdata(:,1)>0),data(sortidx(sortdata(:,1)>0),:));
hold on;
plot([0 0 NaN 2 2 ],[1 sum(sortdata(:,1)>0) NaN 1 sum(sortdata(:,1)>0)],'k:');
h(4) = axes('Position',axpt(10,1,9,1,axpt(2,1,2,1,[],[0.1 0.05])));
imagesc(sortdata(sortdata(:,1)>0,1));
colormap(h(4),ct);
h(6) = axes('Position',axpt(10,1,10,1,axpt(2,1,2,1,[],[0.1 0.05])));
imagesc(sortdata(sortdata(:,1)>0,3));
colormap(h(6),ct2);

set(h(4),'CLim',[0 3])
set(h,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35);
set(h([1,3]),'CLim',[-2 6],'XLim',[-0.5 2.5],'XTick',0:1:2)
set(h(1:2),'YTick',0:500:2500);
set(h(3:4),'YTick',0:200:800);
set(h([2,4,5,6]),'YTickLabel',[],'XTick',[]);
cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig3');
print(fHandle,'-depsc','-painters','psth_gratings.ai');

%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.15 20*0.18]);
hc = colorbar('North');
set(gca,'CLim',[-2 6]);
set(hc,'Box','off','TickDir','out','FontSize',6,'XTick',-2:2:6);
axis off
cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig3');
print(fHandle,'-depsc','-painters','colorbar.ai')

%%
data(isinf(data)) = NaN;
ct = [0.6 0.6 0.6;cbrewer('qual','Dark2',3)];
close all
clear h
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 7]);
h(1) = axes('Position',axpt(1,2,1,1));
hold on;
for iClst = 1:nclst+1
    in = sortdata(:,1)==iClst-1 & sortdata(:,2)>0;
    m = nanmean(data(sortidx(in),:));
    s = nanstd(data(sortidx(in),:))/sqrt(sum(in));
    fill([fr_gratings.psth.time flip(fr_gratings.psth.time)],[m+s flip(m-s)],...
        ct(iClst,:),'EdgeColor','none');
    plot(fr_gratings.psth.time,m,'Color',ct(iClst,:));
end
plot([0 0 NaN 2 2],[-1 5 NaN -1 5],'k:');
alpha(0.2);

h(2) = axes('Position',axpt(1,2,1,2));
hold on;
for iClst = 1:nclst+1
    in = sortdata(:,1)==iClst-1 & sortdata(:,2)<0;
    m = nanmean(data(sortidx(in),:));
    s = nanstd(data(sortidx(in),:))/sqrt(sum(in));
    fill([fr_gratings.psth.time flip(fr_gratings.psth.time)],[m+s flip(m-s)],...
        ct(iClst,:),'EdgeColor','none');
    plot(fr_gratings.psth.time,m,'Color',ct(iClst,:));
end
plot([0 0 NaN 2 2],[-1 5 NaN -1 5],'k:');
set(h,'XLim',[-0.5 2.5],'YLim',[-1 5]);
alpha(0.2);

% xlabel('

%%
close all

sig = cell2mat(cellfun(@(x,y) x&y<0.01,sigrf,p_responsive,'UniformOutput',false));
clusterp = cell2mat(cluster);
clusterp(isnan(clusterp)) = 0;
areap = cat(1,area{:});
nCell = 100;
nIter = 100;

areaLabel = {'V1';'LM';'RL';'AL';'PM';'AM'};
areaList = {'VISp';'VISl';'VISrl';'VISal';'VISpm';'VISam'};
nArea = length(areaList);
datasub = cell(nclst+1,nArea);
ct = [0.6 0.6 0.6;cbrewer('qual','Dark2',3)];

intest = cellfun(@(x,z) x&ismember(z,areaList),sigrf,area,'UniformOutput',false);

group1 = cellfun(@(x,y) x(y),cluster,intest,'UniformOutput',false);
group2 = cellfun(@(x,y) x(y),area,intest,'UniformOutput',false);
% group3 = cellfun(@(x,y) x(y)>0,responsiveness,intest,'UniformOutput',false);
group = {cell2mat(group1),categorical(cat(1,group2{:}))};

for iD = 1:4
    switch iD
        case 1
            data = cell2mat(responsiveness);
            binrange = -2:0.1:5;
            ylimit = [-0.5 1.5];
            xticks = -2:2:4;
            datalist = 'Responsiveness';
        case 2
            data = abs(cell2mat(selectivity));
            binrange = 0:0.1:5;
            ylimit = [0 1.2];
            xticks = 0:1:5;
            datalist = 'Selectivity';
        case 3            
            data = log10(cell2mat(modindex));
            binrange = -3:0.1:1.5;
            ylimit = [-1 1];
            xticks = -3:1:1;
            datalist = 'log10(MI)';
        case 4
            data = cell2mat(latency);
            binrange = 0:0.01:1.5;
            xticks = 0:0.5:1.5;
            ylimit = [0 0.4];
            datalist = 'Latency (s)';
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
    text(1,diff(ylimit)*0.95+ylimit(1),['P_c_l_u_s_t_e_r = ',num2str(round(panova(1)*1000)/1000)],'FontSize',6);
    text(1,diff(ylimit)*0.87+ylimit(1),['P_a_r_e_a = ',num2str(round(panova(2)*1000)/1000)],'FontSize',6);
    text(1,diff(ylimit)*0.79+ylimit(1),['P_c_X_a = ',num2str(round(panova(3)*1000)/1000)],'FontSize',6);
    
    cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig3');
    print(fHandle,'-depsc','-painters',['dg_',datalist,'.ai']);
end

%% 

group = {cell2mat(group1),categorical(cat(1,group2{:}))};

r = cell2mat(cellfun(@(x,y) x(y),responsiveness,in,'UniformOutput',false));
s = cell2mat(cellfun(@(x,y) x(y),selectivity,in,'UniformOutput',false));
m = cell2mat(cellfun(@(x,y) log10(x(y)),modindex,in,'UniformOutput',false));
l = cell2mat(cellfun(@(x,y) x(y),latency,in,'UniformOutput',false));

anovan(r,group,'model','interaction','display','off');
anovan(s,group,'model','interaction');
anovan(m,group,'model','interaction');
anovan(l,group,'model','interaction');





