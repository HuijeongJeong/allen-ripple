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

Tfst = readtable('D:\OneDrive\1.allen-andermann\time_to_first_spike.csv');
Tts = readtable('D:\OneDrive\1.allen-andermann\unit_table_autocorr.csv');
[pval_rf,ctype,cluster,area,tfs,rfsize,timescale,modidx,fr_dg] = deal(cell(nS,1));

for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    
    invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis));
    inthal = ismember(T.unit_id,tag.info.unit_id(ismember(tag.info.structure,{'LGd';'LGv';'LP'})));
    
    [in,idx] = ismember(T.unit_id,unitid);
    clstidx = NaN(length(T.unit_id),1);
    clstidx(in) = cluster_idx.vis{nclst-1}(idx(in));
    
    [~,idx] = ismember(T.unit_id,tag.info.unit_id);
    celltype = [tag.celltype.rs(idx), tag.celltype.fs(idx)];
    
    in = invis | inthal;
  
    ctype{iS} = celltype(in,:);
    cluster{iS} = clstidx(in);
    area{iS} = T.ecephys_structure_acronym(in);
    pval_rf{iS} = T.p_value_rf(in);
    fr_dg{iS} = T.firing_rate_dg(in);
        
    [~,idx] = ismember(T.unit_id(in),Tfst.ecephys_unit_id);
    tfs{iS} = Tfst.time_to_first_spike_fl(idx);
    rfsize{iS} = T.area_rf(in);
    modidx{iS} = T.mod_idx_dg(in);
        
    timescale{iS} = NaN(sum(in),1);
    [inn,idx] = ismember(T.unit_id(in),Tts.unit_id);
    timescale{iS}(inn) = Tts.timescale_ac(idx(inn));
    intimescale = Tts.err_ac(in)<20 & Tts.spike_count_ac(in)>50;
    timescale{iS}(~intimescale) = nan;
end

ct = [0 0 0; cbrewer('qual','Dark2',3)];
binrange = {[0:50:2500];[30:3:100];[-2:0.1:1];[0:3:150]};
measureList = {'Receptive field area (deg^2)','Time to first spike (ms)',...
    'log_1_0(MI)','Response decay timescale (ms)'};
iCT = 1;

areaList = {'VISp';'VISl';'VISrl';'VISal';'VISpm';'VISam'};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 15 7]);
for iClst = 1:4
    if iClst==1
        in = cellfun(@(x,y,z,w,v) isnan(x) & y(:,iCT) & z<0.01 & w<2500 & v<0.1,...
            cluster,ctype,pval_rf,rfsize,tfs,'UniformOutput',false);
    else
        in = cellfun(@(x,y,z,w,v) x==iClst-1 & y(:,iCT) & z<0.01 & w<2500 & v<0.1,...
            cluster,ctype,pval_rf,rfsize,tfs,'UniformOutput',false);
    end
    data = cell2mat(cellfun(@(x,y,z,w,v) [x(v),y(v),z(v),w(v)],rfsize,tfs,modidx,timescale,in,'UniformOutput',false));

    areadata = cellfun(@(x,y) x(y),area,in,'UniformOutput',false);
    areadata = cat(1,areadata{:});
    for iMs = 1:4
        if iClst==1
            h(iMs) = axes('Position',axpt(4,2,iMs,1,[],[0.1 0.15])); hold on;
            h(iMs+4) = axes('Position',axpt(4,2,iMs,2,[],[0.1 0.15])); hold on;
        end
        if iMs==3
            bincount = histc(log10(data(:,iMs)),binrange{iMs});
        elseif iMs==2
            bincount = histc((data(:,iMs))*1000,binrange{iMs});
        else
            bincount = histc(data(:,iMs),binrange{iMs});
        end
        plot(binrange{iMs},cumsum(bincount)./sum(bincount),...
            'Color',ct(iClst,:),'Parent',h(iMs));
        [m,s] = deal(NaN(length(areaList),1));
        for iA = 1:length(areaList)
            sdata = data(strcmp(areadata,areaList{iA}),iMs);
            if iMs==3
                sdata = log10(sdata);
            elseif iMs==2
                sdata = sdata*1000;
            end
            m(iA) = nanmean(sdata);
            s(iA) = nanstd(sdata)/sqrt(size(sdata,1));
        end
        errorbar(1:length(areaList),m,s,'Color',ct(iClst,:),'Parent',h(iMs+4),'CapSize',3)
        %                     errorbar(iClst+(iA-1)*6,nanmean(sdata),nanstd(sdata)/sqrt(size(sdata,1)),'Color',ct(iClst,:),'Parent',h(iMs+4));

        if iClst==4
            set(h(iMs),'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,'XLim',[min(binrange{iMs}),max(binrange{iMs})])
            set(h(iMs+4),'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35,...
                'XLim',[0 length(areaList)+1],'XTick',1:6,'XTickLabel',areaList,...
                'XTickLabelRotation',45);
            xlabel(h(iMs),measureList{iMs});
            ylabel(h(iMs+4),measureList{iMs});
            ylabel(h(iMs),'cumulative fraction')
        end
    end
end
cd('D:\OneDrive - UCSF\figures\allen\cluster_characteristic')
print(fHandle,'-dtiff','-r600','cluster_hierarchy_characteristics')
