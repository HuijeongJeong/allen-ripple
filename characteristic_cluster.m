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

stimList = {'dot_motion';'flashes';'gabors'};
[firingrate,laminar_loc,sparseness,fanofactor,else_dg,timetopeak,adapt_r,adapt_p,...
    area_rf,pval_rf,behavior_corr,behavior_corrp,ctype,cluster,area,history,...
    behavior_selec,behavior_selecp] = deal(cell(nS,1));
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_stimulus_evoked.mat'],'fr_behavior','fr_movie');
    load([sdir(sessionList(iS)),'_Events.mat'],'drifting_gratings','stimulus_presentation');
    load([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_win');
    
    invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis));
    inthal = ismember(T.unit_id,tag.info.unit_id(ismember(tag.info.structure,{'LGd';'LGv';'LP'})));
    
    [in,idx] = ismember(T.unit_id,unitid);
    clstidx = NaN(length(T.unit_id),1);
    clstidx(in) = cluster_idx.vis{nclst-1}(idx(in));
    
    [~,idx] = ismember(T.unit_id,tag.info.unit_id);
    celltype = [tag.celltype.rs(idx), tag.celltype.fs(idx)];
    laminar_loc{iS} = tag.info.distance_from_L4(idx);
    
    in = invis | inthal;
  
    laminar_loc{iS} = laminar_loc{iS}(in);
    ctype{iS} = celltype(in,:);
    cluster{iS} = clstidx(in);
    area{iS} = T.ecephys_structure_acronym(in);
    
    fr_spontaneous = cellfun(@(x) sum(x>=spontaneous_win(1) & x<=spontaneous_win(2))/...
        diff(spontaneous_win),T.spike_time);
    firingrate{iS} = [T.firing_rate_dg(in), T.firing_rate_dm(in),...
        T.firing_rate_fl(in), T.firing_rate_rf(in), fr_movie.psth.mean_total(in),fr_spontaneous(in)];
    sparseness{iS} = [T.lifetime_sparseness_dg(in), T.lifetime_sparseness_dm(in),...
        T.lifetime_sparseness_fl(in), T.lifetime_sparseness_rf(in)];
    fanofactor{iS} = [T.fano_dg(in), T.fano_dm(in), T.fano_fl(in), T.fano_rf(in)];
    timetopeak{iS} = T.time_to_peak_dm(in);
    else_dg{iS} = [T.mod_idx_dg(in), T.f1_f0_dg(in), T.c50_dg(in)];
    
    [adapt_r{iS},adapt_p{iS}] = deal(NaN(sum(in),4));
    history{iS} = NaN(sum(in),2);
    [adapt_r{iS}(:,1),adapt_p{iS}(:,1)] = adaptationidex(T.spike_time(in),drifting_gratings.window,...
        drifting_gratings.blockInd(:,1),drifting_gratings.stimulus_condition_id);  
    history{iS}(:,1) = historyidx(T.spike_time(in),drifting_gratings.window,...
        drifting_gratings.blockInd(:,1),drifting_gratings.stimulus_condition_id);  
    
    for iStim = 1:3
        instim = strcmp(stimulus_presentation.stimulus_name,stimList{iStim});
        window = [stimulus_presentation.start_time(instim), stimulus_presentation.stop_time(instim)];
        condition_id = stimulus_presentation.stimulus_condition_id(instim);
        [adapt_r{iS}(:,iStim+1),adapt_p{iS}(:,iStim+1)] = adaptationidex(T.spike_time(in),...
            window,ones(size(window,1),1),condition_id);
        if iStim==2
            history{iS}(:,iStim) = historyidx(T.spike_time(in),window,ones(size(window,1),1),condition_id);
        end
    end
    
    area_rf{iS} = T.area_rf(in);
    pval_rf{iS} = T.p_value_rf(in);
    
    behavior_corr{iS} = [fr_behavior.run.corr.r(in), fr_behavior.pupil.corr.r(in)];
    behavior_corrp{iS} = [fr_behavior.run.corr.pval(in), fr_behavior.pupil.corr.pval(in)];
    
    behavior_selec{iS} = [fr_behavior.run.selec.idx(in), fr_behavior.pupil.selec.idx(in)];
    behavior_selecp{iS} = [fr_behavior.run.selec.pval(in), fr_behavior.pupil.selec.pval(in)];
end

ct = [0.8 0.8 0.8; cbrewer('qual','Pastel2',3); 0,0,0; cbrewer('qual','Dark2',3)];
visList = {'VISp';'VISl';'VISal';'VISrl';'VISam';'VISpm'};
celltype = {'RS';'FS'};
nVis = length(visList);

%% behavior plot
behList = {'running';'pupil'};

for ii = 1:2
    if ii==1
        p = behavior_corrp;
        idx = behavior_corr;
        datatype = 'corr';
        binrange = -0.6:0.05:0.6;
    else
        p = behavior_selecp;
        idx = behavior_selec;
        datatype = 'selec';
        binrange = -1.5:0.05:1.5;
    end
    
    fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 6]);
    for iCT = 1:2
        h(1+(2-1)*iCT) = axes('Position',axpt(2,2,1,iCT,[],[0.2 0.2])); hold on;
        h(2+(2-1)*iCT) = axes('Position',axpt(2,2,2,iCT,[],[0.2 0.2])); hold on;
        for iClst = 1:4
           for iVis = 1:nVis
                if iClst==1
                    in = cellfun(@(x,y,z,w) isnan(x) & y(:,iCT) & strcmp(z,visList{iVis}) & w<0.01,...
                        cluster,ctype,area,pval_rf,'UniformOutput',false);
                else
                    in = cellfun(@(x,y,z,w) x==iClst-1 & y(:,iCT) & strcmp(z,visList{iVis}) & w<0.01,...
                        cluster,ctype,area,pval_rf,'UniformOutput',false);
                end
                
                data_p = cell2mat(cellfun(@(x,y) x(y,:),p,in,'UniformOutput',false));
                data_r = cell2mat(cellfun(@(x,y) x(y,:),idx,in,'UniformOutput',false));
                
                for ibeh = 1:2
                    bincount = histc(data_r(:,ibeh),binrange);
                    data{ibeh,iVis,iClst} = cumsum(bincount)/sum(bincount);
                    plot(binrange,data{ibeh,iVis,iClst},'Color',ct(iClst,:),'Parent',h(ibeh+(2-1)*iCT));
                    fon(iVis,ibeh,iClst) = nanmean(data_p(:,ibeh)<0.05);
                end
            end
        end
        
        for iClst = 1:4
            for ibeh = 1:2
                plot(binrange,nanmean(cell2mat(data(ibeh,:,iClst)),2),'Color',ct(iClst+4,:),'Parent',h(ibeh+(2-1)*iCT));
                text((binrange(end)-binrange(1))*0.1+binrange(1),...
                    0.95-(iClst-1)*0.07,[num2str(round(nanmean(fon(:,ibeh,iClst))*1000)/10),'%'],...
                    'Color',ct(iClst+4,:),'Parent',h(ibeh+(2-1)*iCT),'FontSize',5);
                if iClst==4
                    plot([0 0],[0 1],'k:','Parent',h(ibeh+(2-1)*iCT));
                    title(h(ibeh+(2-1)*iCT),[celltype{iCT},'-',behList{ibeh}]);
                    xlabel(h(ibeh+(2-1)*iCT),'r','FontSize',5);
                    ylabel(h(ibeh+(2-1)*iCT),'cumulative fraction','FontSize',5);
                    set(h(ibeh+(2-1)*iCT),'Box','off','TickDir','out',...
                        'FontSize',5,'XLim',[binrange(1) binrange(end)]);
                end
            end
        end
    end
    cd('D:\OneDrive - University of California, San Francisco\figures\allen\cluster_characteristic');
    print(fHandle,'-dtiff','-r600',['behavior_response',datatype,'.tif']);
end

%% receptive field
clear h data fon
binrange = 0:100:4500;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 6]);
for iCT = 1:2
    h(iCT) = axes('Position',axpt(10,2,2:10,iCT,[],[0.2 0.2])); hold on;
    for iClst = 1:4
        for iVis = 1:nVis
            if iClst==1
                in = cellfun(@(x,y,z) isnan(x) & y(:,iCT) & strcmp(z,visList{iVis}),...
                    cluster,ctype,area,'UniformOutput',false);
            else
                in = cellfun(@(x,y,z) x==iClst-1 & y(:,iCT) & strcmp(z,visList{iVis}),...
                    cluster,ctype,area,'UniformOutput',false);
            end
            data_p = cell2mat(cellfun(@(x,y) x(y),pval_rf,in,'UniformOutput',false));
            data_r = cell2mat(cellfun(@(x,y) x(y),area_rf,in,'UniformOutput',false));
            
            bincount = histc(data_r(data_p<0.01),binrange);
            data{iVis,iClst} = cumsum(bincount)/sum(bincount);
            plot(binrange,data{iVis,iClst},'Color',ct(iClst,:),'Parent',h(iCT));
            fon(iVis,iClst) = nanmean(data_p<0.01);
        end
    end
    
    for iClst = 1:4
        plot(binrange,nanmean(cell2mat(data(:,iClst)'),2),'Color',ct(iClst+4,:),'Parent',h(iCT));
        text(4500*0.8,0.4-(iClst-1)*0.07,[num2str(round(nanmean(fon(:,iClst))*1000)/10),'%'],...
            'Color',ct(iClst+4,:),'Parent',h(iCT),'FontSize',5);
        if iClst==4
            title(h(iCT),[celltype{iCT}]);
            xlabel(h(iCT),'area (receptive field)','FontSize',5);
            ylabel(h(iCT),'cumulative fraction','FontSize',5);
            set(h(iCT),'Box','off','TickDir','out','FontSize',5,'XLim',[0 4500]);
        end
    end
end
cd('D:\OneDrive - University of California, San Francisco\figures\allen\cluster_characteristic');
print(fHandle,'-dtiff','-r600','receptive_field.tif');

%% sparsness
clear h data fon
stimList = {'drifting grating';'dot motion';'flash';'gabors'};
binrange = 0:0.01:1;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 15 6]);
for iCT = 1:2
   h(1+(2-1)*iCT) = axes('Position',axpt(4,2,1,iCT,[],[0.1 0.2])); hold on;
   h(2+(2-1)*iCT) = axes('Position',axpt(4,2,2,iCT,[],[0.1 0.2])); hold on;
   h(3+(2-1)*iCT) = axes('Position',axpt(4,2,3,iCT,[],[0.1 0.2])); hold on;
   h(4+(2-1)*iCT) = axes('Position',axpt(4,2,4,iCT,[],[0.1 0.2])); hold on;
   for iClst = 1:4
       for iVis = 5
           if iClst==1
               in = cellfun(@(x,y,z) isnan(x) & y(:,iCT) & strcmp(z,visList{iVis}),...
                   cluster,ctype,area,'UniformOutput',false);
           else
               in = cellfun(@(x,y,z) x==iClst-1 & y(:,iCT) & strcmp(z,visList{iVis}),...
                   cluster,ctype,area,'UniformOutput',false);
           end
           data_r = cell2mat(cellfun(@(x,y) x(y,:),sparseness,in,'UniformOutput',false));
           data_p = cell2mat(cellfun(@(x,y) x(y),pval_rf,in,'UniformOutput',false));
           for istim = 1:4
              bincount = histc(data_r(data_p<0.01,istim),binrange);
              data{istim,iVis,iClst} = cumsum(bincount)/sum(bincount);
              plot(binrange,data{istim,iVis,iClst},'Color',ct(iClst,:),'Parent',h(istim+(2-1)*iCT));
           end
       end
   end
   
   for iClst = 1:4
       for istim = 1:4
           plot(binrange,nanmean(cell2mat(data(istim,:,iClst)),2),'Color',ct(iClst+4,:),'Parent',h(istim+(2-1)*iCT));
           if iClst==4
              title(h(istim+(2-1)*iCT),[celltype{iCT},'-',stimList{istim}]);
              xlabel(h(istim+(2-1)*iCT),'lifetime sparseness','FontSize',5);
              ylabel(h(istim+(2-1)*iCT),'cumulative fraction','FontSize',5);
              set(h(istim+(2-1)*iCT),'Box','off','TickDir','out','FontSize',5);
           end
       end
   end
end
cd('D:\OneDrive - University of California, San Francisco\figures\allen\cluster_characteristic');
print(fHandle,'-dtiff','-r600','sparseness.tif');

%% fano factor
clear h data fon
stimList = {'drifting grating';'dot motion';'flash';'gabors'};
binrange = 0:0.5:30;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 15 6]);
for iCT = 1:2
   h(1+(2-1)*iCT) = axes('Position',axpt(4,2,1,iCT,[],[0.1 0.2])); hold on;
   h(2+(2-1)*iCT) = axes('Position',axpt(4,2,2,iCT,[],[0.1 0.2])); hold on;
   h(3+(2-1)*iCT) = axes('Position',axpt(4,2,3,iCT,[],[0.1 0.2])); hold on;
   h(4+(2-1)*iCT) = axes('Position',axpt(4,2,4,iCT,[],[0.1 0.2])); hold on;
   for iClst = 1:4
       for iVis = nVis
           if iClst==1
               in = cellfun(@(x,y,z) isnan(x) & y(:,iCT) & strcmp(z,visList{iVis}),...
                   cluster,ctype,area,'UniformOutput',false);
           else
               in = cellfun(@(x,y,z) x==iClst-1 & y(:,iCT) & strcmp(z,visList{iVis}),...
                   cluster,ctype,area,'UniformOutput',false);
           end
           data_r = cell2mat(cellfun(@(x,y) x(y,:),fanofactor,in,'UniformOutput',false));
           data_p = cell2mat(cellfun(@(x,y) x(y),pval_rf,in,'UniformOutput',false));
           for istim = 1:4
              bincount = histc(data_r(data_p<0.01,istim),binrange);
              data{istim,iVis,iClst} = cumsum(bincount)/sum(bincount);
              plot(binrange,data{istim,iVis,iClst},'Color',ct(iClst,:),'Parent',h(istim+(2-1)*iCT));
           end
       end
   end
   
   for iClst = 1:4
       for istim = 1:4
           plot(binrange,nanmean(cell2mat(data(istim,:,iClst)),2),'Color',ct(iClst+4,:),'Parent',h(istim+(2-1)*iCT));
           if iClst==4
              title(h(istim+(2-1)*iCT),[celltype{iCT},'-',stimList{istim}]);
              xlabel(h(istim+(2-1)*iCT),'fano factor','FontSize',5);
              ylabel(h(istim+(2-1)*iCT),'cumulative fraction','FontSize',5);
              set(h(istim+(2-1)*iCT),'Box','off','TickDir','out','FontSize',5);
           end
       end
   end
end
cd('D:\OneDrive - University of California, San Francisco\figures\allen\cluster_characteristic');
print(fHandle,'-dtiff','-r600','fano_factor.tif');

%% other-drifting gratings
clear h data fon
binrange = {0:0.1:10;0:0.01:2;0:0.005:1};
idxList = {'mod index';'f1/f0';'c50'};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 11 6]);
for iCT = 1:2
   h(1+(2-1)*iCT) = axes('Position',axpt(3,2,1,iCT,[],[0.1 0.2])); hold on;
   h(2+(2-1)*iCT) = axes('Position',axpt(3,2,2,iCT,[],[0.1 0.2])); hold on;
   h(3+(2-1)*iCT) = axes('Position',axpt(3,2,3,iCT,[],[0.1 0.2])); hold on;
   for iClst = 1:4
       for iVis = 3
           if iClst==1
               in = cellfun(@(x,y,z) isnan(x) & y(:,iCT) & strcmp(z,visList{iVis}),...
                   cluster,ctype,area,'UniformOutput',false);
           else
               in = cellfun(@(x,y,z) x==iClst-1 & y(:,iCT) & strcmp(z,visList{iVis}),...
                   cluster,ctype,area,'UniformOutput',false);
           end
           data_p = cell2mat(cellfun(@(x,y) x(y),pval_rf,in,'UniformOutput',false));
           data_r = cell2mat(cellfun(@(x,y) x(y,:),else_dg,in,'UniformOutput',false));
           
           for ibeh = 1:3
              bincount = histc(data_r(data_p<0.01,ibeh),binrange{ibeh});
              data{ibeh,iVis,iClst} = cumsum(bincount)/sum(bincount);
              plot(binrange{ibeh},data{ibeh,iVis,iClst},'Color',ct(iClst,:),'Parent',h(ibeh+(2-1)*iCT));
           end
       end
   end
   
   for iClst = 1:4
       for ibeh = 1:3
           plot(binrange{ibeh},nanmean(cell2mat(data(ibeh,:,iClst)),2),'Color',ct(iClst+4,:),'Parent',h(ibeh+(2-1)*iCT));
           if iClst==4
               if ibeh==2
                  plot([1 1],[0 1],'k:','Parent',h(ibeh+(2-1)*iCT));
               end
              title(h(ibeh+(2-1)*iCT),celltype{iCT});
              xlabel(h(ibeh+(2-1)*iCT),idxList{ibeh},'FontSize',5);
              ylabel(h(ibeh+(2-1)*iCT),'cumulative fraction','FontSize',5);
              set(h(ibeh+(2-1)*iCT),'Box','off','TickDir','out','FontSize',5);
           end
       end
   end
end
cd('D:\OneDrive - University of California, San Francisco\figures\allen\cluster_characteristic');
print(fHandle,'-dtiff','-r600','dg_response.tif');

%% firing rate
clear h data data2 fon

stimList = {'drifting grating';'dot motion';'flash';'gabors';'movie';'spontaneous'};
binrange = 0:0.5:50;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 6]);
for iCT = 1:2
   h(1+(2-1)*iCT) = axes('Position',axpt(2,2,1,iCT,[],[0.2 0.2])); hold on;
   h(2+(2-1)*iCT) = axes('Position',axpt(2,2,2,iCT,[],[0.2 0.2])); hold on;
   
   for iClst = 1:4
       for iVis = 1:nVis
           if iClst==1
               in = cellfun(@(x,y,z) isnan(x) & y(:,iCT) & strcmp(z,visList{iVis}),...
                   cluster,ctype,area,'UniformOutput',false);
           else
               in = cellfun(@(x,y,z) x==iClst-1 & y(:,iCT) & strcmp(z,visList{iVis}),...
                   cluster,ctype,area,'UniformOutput',false);
           end
           data_r = cell2mat(cellfun(@(x,y) x(y,:),firingrate,in,'UniformOutput',false));
           data_p = cell2mat(cellfun(@(x,y) x(y),pval_rf,in,'UniformOutput',false));
           bincount = histc(data_r(data_p<0.01,6),binrange);
           data{iVis,iClst} = cumsum(bincount)/sum(bincount);
           plot(binrange,data{iVis,iClst},'Color',ct(iClst,:),'Parent',h(1+(2-1)*iCT));
           
           data2{iClst,iCT}(iVis,:) = nanmean(data_r(data_p<0.05,:));
       end
   end
   
   for iClst = 1:4
       plot(binrange,nanmean(cell2mat(data(:,iClst)'),2),'Color',ct(iClst+4,:),'Parent',h(1+(2-1)*iCT));
       plot(1:6,data2{iClst,iCT}','Color',ct(iClst,:),'Parent',h(2+(2-1)*iCT));
       plot(1:6,nanmean(data2{iClst,iCT}),'Color',ct(iClst+4,:),'Parent',h(2+(2-1)*iCT));
       if iClst==4
           title(h(1+(2-1)*iCT),[celltype{iCT},'-',stimList{6}]);
           xlabel(h([1:2]+(2-1)*iCT),'firing rate','FontSize',5);
           ylabel(h(1+(2-1)*iCT),'cumulative fraction','FontSize',5);
           set(h([1:2]+(2-1)*iCT),'Box','off','TickDir','out','FontSize',5);
           set(h(2+(2-1)*iCT),'XTick',1:6,'XTickLabel',stimList,'XTickLabelRotation',45,'XLim',[0 7]);
       end
   end
end
cd('D:\OneDrive - University of California, San Francisco\figures\allen\cluster_characteristic');
print(fHandle,'-dtiff','-r600','firing_rate.tif');

%% layer distribution
clear h data fon
binrange = -400:10:800;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 6]);
for iCT = 1:2
    h(iCT) = axes('Position',axpt(10,2,2:10,iCT,[],[0.2 0.2])); hold on;
    for iClst = 1:4
        for iVis = 1:nVis
            if iClst==1
                in = cellfun(@(x,y,z) isnan(x) & y(:,iCT) & strcmp(z,visList{iVis}),...
                    cluster,ctype,area,'UniformOutput',false);
            else
                in = cellfun(@(x,y,z) x==iClst-1 & y(:,iCT) & strcmp(z,visList{iVis}),...
                    cluster,ctype,area,'UniformOutput',false);
            end
            data_p = cell2mat(cellfun(@(x,y) x(y),pval_rf,in,'UniformOutput',false));
            data_r = cell2mat(cellfun(@(x,y) x(y),laminar_loc,in,'UniformOutput',false));
            
            bincount = histc(data_r(data_p<0.01),binrange);
            data{iVis,iClst} = cumsum(bincount)/sum(bincount);
            plot(binrange,data{iVis,iClst},'Color',ct(iClst,:),'Parent',h(iCT));
            fon(iVis,iClst) = nanmean(data_p<0.01);
        end
    end
    
    for iClst = 1:4
        plot(binrange,nanmean(cell2mat(data(:,iClst)'),2),'Color',ct(iClst+4,:),'Parent',h(iCT));
        
        if iClst==4
            title(h(iCT),[celltype{iCT}]);
            xlabel(h(iCT),'distance from the first sink (um)','FontSize',5);
            ylabel(h(iCT),'cumulative fraction','FontSize',5);
            set(h(iCT),'Box','off','TickDir','out','FontSize',5,'XLim',[-400 800]);
        end
    end
end
cd('D:\OneDrive - University of California, San Francisco\figures\allen\cluster_characteristic');
print(fHandle,'-dtiff','-r600','laminar_location.tif');

%% adaptation index
stimList = {'drifting grating';'dot motion';'flash';'gabors'};
binrange = -1:0.05:1;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 15 6]);
for iCT = 1:2
   h(1+(2-1)*iCT) = axes('Position',axpt(4,2,1,iCT,[],[0.1 0.2])); hold on;
   h(2+(2-1)*iCT) = axes('Position',axpt(4,2,2,iCT,[],[0.1 0.2])); hold on;
   h(3+(2-1)*iCT) = axes('Position',axpt(4,2,3,iCT,[],[0.1 0.2])); hold on;
   h(4+(2-1)*iCT) = axes('Position',axpt(4,2,4,iCT,[],[0.1 0.2])); hold on;
   for iClst = 1:4
       for iVis = 1:nVis
           if iClst==1
               in = cellfun(@(x,y,z) isnan(x) & y(:,iCT) & strcmp(z,visList{iVis}),...
                   cluster,ctype,area,'UniformOutput',false);
           else
               in = cellfun(@(x,y,z) x==iClst-1 & y(:,iCT) & strcmp(z,visList{iVis}),...
                   cluster,ctype,area,'UniformOutput',false);
           end
           data_r = cell2mat(cellfun(@(x,y) x(y,:),adapt_r,in,'UniformOutput',false));
           data_p = cell2mat(cellfun(@(x,y) x(y),pval_rf,in,'UniformOutput',false));
           for istim = 1:4
              bincount = histc(data_r(data_p<0.01,istim),binrange);
              data{istim,iVis,iClst} = cumsum(bincount)/sum(bincount);
              plot(binrange,data{istim,iVis,iClst},'Color',ct(iClst,:),'Parent',h(istim+(2-1)*iCT));
           end
       end
   end
   
   for iClst = 1:4
       for istim = 1:4
           plot(binrange,nanmean(cell2mat(data(istim,:,iClst)),2),'Color',ct(iClst+4,:),'Parent',h(istim+(2-1)*iCT));
           if iClst==4
              title(h(istim+(2-1)*iCT),[celltype{iCT},'-',stimList{istim}]);
              xlabel(h(istim+(2-1)*iCT),'adaptation index','FontSize',5);
              ylabel(h(istim+(2-1)*iCT),'cumulative fraction','FontSize',5);
              set(h(istim+(2-1)*iCT),'Box','off','TickDir','out','FontSize',5);
           end
       end
   end
end
cd('D:\OneDrive - University of California, San Francisco\figures\allen\cluster_characteristic');
print(fHandle,'-dtiff','-r600','adaptation_index.tif');

%% history index
behList = {'drifting grating';'flash'};
binrange = -1:0.05:1;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 6]);
for iCT = 1:2
   h(1+(2-1)*iCT) = axes('Position',axpt(2,2,1,iCT,[],[0.2 0.2])); hold on;
   h(2+(2-1)*iCT) = axes('Position',axpt(2,2,2,iCT,[],[0.2 0.2])); hold on;
   for iClst = 1:4
       for iVis = 1:nVis
           if iClst==1
               in = cellfun(@(x,y,z,w) isnan(x) & y(:,iCT) & strcmp(z,visList{iVis}) & w<0.01,...
                   cluster,ctype,area,pval_rf,'UniformOutput',false);
           else
               in = cellfun(@(x,y,z,w) x==iClst-1 & y(:,iCT) & strcmp(z,visList{iVis}) & w<0.01,...
                   cluster,ctype,area,pval_rf,'UniformOutput',false);
           end
           data_r = cell2mat(cellfun(@(x,y) x(y,:),history,in,'UniformOutput',false));
           
           for ibeh = 1:2
              bincount = histc(data_r(:,ibeh),binrange);
              data{ibeh,iVis,iClst} = cumsum(bincount)/sum(bincount);
              plot(binrange,data{ibeh,iVis,iClst},'Color',ct(iClst,:),'Parent',h(ibeh+(2-1)*iCT));
           end
       end
   end
   
   for iClst = 1:4
       for ibeh = 1:2
           plot(binrange,nanmean(cell2mat(data(ibeh,:,iClst)),2),'Color',ct(iClst+4,:),'Parent',h(ibeh+(2-1)*iCT));
           if iClst==4
              plot([0 0],[0 1],'k:','Parent',h(ibeh+(2-1)*iCT));
              title(h(ibeh+(2-1)*iCT),[celltype{iCT},'-',behList{ibeh}]);
              xlabel(h(ibeh+(2-1)*iCT),'history index','FontSize',5);
              ylabel(h(ibeh+(2-1)*iCT),'cumulative fraction','FontSize',5);
              set(h(ibeh+(2-1)*iCT),'Box','off','TickDir','out','FontSize',5);
           end
       end
   end
end
cd('D:\OneDrive - University of California, San Francisco\figures\allen\cluster_characteristic');
print(fHandle,'-dtiff','-r600','history_dependency.tif');


function [radap,padap] = adaptationidex(spiketime,window,block,stimulustype)
fr = cellfun(@(x) cellfun(@(y) sum(x>=y(1) & x<=y(2))/diff(y),...
    mat2cell(window,ones(size(window,1),1),2)),spiketime,'UniformOutput',false);
stimuluslist = unique(stimulustype);
avefr = NaN(length(spiketime),length(stimuluslist));
for iS = 1:length(stimuluslist)
   avefr(:,iS) = cellfun(@(x) nanmean(x(stimulustype==stimuluslist(iS))),fr);
end
[~,maxidx] = max(avefr,[],2);
[radap,padap] = deal(NaN(length(spiketime),size(block,2)));
for iB = 1:size(block,2)
    [radap(:,iB),padap(:,iB)] = cellfun(@(x,y) corr([1:sum(block(:,iB)&stimulustype==stimuluslist(y))]',...
        x(block(:,iB)&stimulustype==stimuluslist(y))),fr,num2cell(maxidx));
end
end

function idx = historyidx(spiketime,window,block,stimulustype)
fr = cellfun(@(x) cellfun(@(y) sum(x>=y(1) & x<=y(2))/diff(y),...
    mat2cell(window,ones(size(window,1),1),2)),spiketime,'UniformOutput',false);
idx = NaN(length(spiketime),size(block,2));
for iB = 1:size(block,2)
    samestim = cellfun(@(x) mean(x([false;stimulustype(1:end-1)==stimulustype(2:end) & block(2:end,iB)])),fr);
    diffstim = cellfun(@(x) mean(x([false;stimulustype(1:end-1)~=stimulustype(2:end) & block(2:end,iB)])),fr);
    idx(:,iB) = (diffstim-samestim)./(samestim+diffstim);
end
end

function sessionDir = sdir(sessionid)

sessionDir = ['D:\OneDrive\1.allen-andermann\session_',num2str(sessionid),...
    '\session_',num2str(sessionid)];
end