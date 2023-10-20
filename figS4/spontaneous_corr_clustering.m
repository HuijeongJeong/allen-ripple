clc; clearvars; close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
session_id = tag.info.session_id(idx);

sessionList = unique(session_id);
nS = length(sessionList);

binsize = 2;
coefbin = -0.3:0.1:0.5;

[inclst,npair] = deal(nan(nS,length(coefbin)+1,3));
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_win','CA1_ripple_classified');
    load([sdir(sessionList(iS)),'_stimulus_evoked.mat'],'fr_gratings');
    load([sdir(sessionList(iS)),'_Events.mat'],'drifting_gratings');
    
    %%
    [in,idx] = ismember(T.unit_id,unit_id.vis);
    inrs = ismember(T.unit_id,tag.info.unit_id(tag.celltype.rs));
    
    spiketimes = T.spike_time(in&inrs);
    cidx = cluster_idx.vis{2}(idx(in&inrs));
    nclst = [sum(cidx==1), sum(cidx==2), sum(cidx==3)];
    if min(nclst)<5
        continue;
    end
    
    %% spontaneous correlation
    
    spkhist = cell2mat(cellfun(@(x) histc(x,spontaneous_win(1):binsize:spontaneous_win(2))',...
        spiketimes,'UniformOutput',false));        
    spkhist = spkhist(:,1:end-1);
    totalripple = sortrows([CA1_ripple_classified.medial(:,[1,3]);...
        CA1_ripple_classified.lateral(:,[1,3]);...
        CA1_ripple_classified.global],1);
    time = spontaneous_win(1)+binsize/2:binsize:spontaneous_win(2)-binsize/2;
    inripple = sum(cell2mat(cellfun(@(y) cellfun(@(x) x-binsize/2<=y(1) & y(2)<x+binsize/2,num2cell(time)),...
        mat2cell(totalripple,ones(size(totalripple,1),1),2),'UniformOutput',false)),1)>0; 
    inripple = inripple';
    spontcoef = corrcoef(spkhist(:,~inripple)');
    
    %% signal/noise correlation
    spkgrating = cell2mat(cellfun(@(y) cellfun(@(x) sum(x>=0 & x<=2),y),...
        fr_gratings.spikeTime(in&inrs),'UniformOutput',false)');
    aversp = nan(sum(in&inrs),8);
    for ic = 1:8
        aversp(:,ic) = mean(spkgrating(drifting_gratings.stimulus_condition_id==2081+ic,:),1);
    end
    signalcoef = corrcoef(aversp');
    noisecoef = corrcoef(spkgrating);
    
    %%
    for id = 1:3
        switch id
            case 1
                data = spontcoef;
            case 2
                data = noisecoef;
            case 3
                data = signalcoef;
        end
        for ic = 1:length(coefbin)+1
            if ic==1
                [i,j] = find(data<=coefbin(1));
            elseif ic==length(coefbin)+1
                [i,j] = find(data>coefbin(ic-1) & data<1);
            else
                [i,j] = find(data>coefbin(ic-1) & data<=coefbin(ic));
            end
            
            %%
            pair = unique(sort([i,j],2),'rows');
            paircidx = cidx(pair);
            if size(paircidx,2)==1
                paircidx = paircidx';
            end
            inclst(iS,ic,id) = sum(paircidx(:,1)==paircidx(:,2));
            npair(iS,ic,id) = size(pair,1);
        end
    end
end

%%
z = 1.96;
close all
plotbin = -0.35:0.1:0.55;
out = sum(isnan(npair(:,:)) | npair(:,:)==0,2)>0;
fraction = inclst./npair;
xlabelList = {'Spontaneous correlation','Noise correlation','Signal correlation'};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 4.5]);
for id = 1:2
    axes('Position',axpt(2,10,id,1:8));
    plot(plotbin,squeeze(fraction(~out,:,id+1)),'Color',[0.6 0.6 0.6],'LineWidth',0.35);
    hold on;
    errorbar(plotbin,mean(squeeze(fraction(~out,:,id+1))),...
        z*std(squeeze(fraction(~out,:,id+1)))/sqrt(sum(~out)),'k','LineWidth',0.5);
    plot([0 0],[0 1],'k:','LineWidth',0.35);
    xlabel(xlabelList{id+1});
    set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',...
        coefbin,'XTickLabel',{'<-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.4','>0.5'},...
        'YTick',0:0.2:1,'XLim',[-0.4 0.6],'XTickLabelRotation',90);
    if id==1
        ylabel('% within cluster');
    else
        set(gca,'YTickLabel',[]);
    end
end
print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\figures\2.allen\revision\FigS_visualresponse\corr_withinclst_visual.ai');

%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4.5 4.5]);
    axes('Position',axpt(10,10,2:10,1:8));
    plot(plotbin,squeeze(fraction(~out,:,1)),'Color',[0.6 0.6 0.6],'LineWidth',0.35);
    hold on;
    errorbar(plotbin,mean(squeeze(fraction(~out,:,1))),...
        std(squeeze(fraction(~out,:,id+1)))/sqrt(sum(~out)),'k','LineWidth',0.5);
    plot([0 0],[0 1],'k:','LineWidth',0.35);
    xlabel(xlabelList{1});
    set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',...
        coefbin,'XTickLabel',{'<-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.4','>0.5'},...
        'YTick',0:0.2:1,'XLim',[-0.4 0.6],'XTickLabelRotation',90);
        ylabel('% within cluster');

print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\figures\2.allen\revision\figS6\corr_withinclst_spontaneous.ai');
%%