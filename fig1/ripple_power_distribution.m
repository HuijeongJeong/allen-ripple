clearvars; clc; close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

refList = {'medial';'lateral';'global_m';'global_l'};
[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
session_id = tag.info.session_id(idx);
sessionList = unique(session_id);
nS = length(sessionList);

binrange = -1:0.1:50;
bincount = cell(2,3);
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_Events.mat'],'probe');
    load([sdir(sessionList(iS)),'_ripples.mat'],'filtered_lfp',...
        'CA1_ripple_classified','spontaneous_win');
    
    for iRp = 1:2 % the probe where signal is from
        if ~iscell(filtered_lfp.lfp{iRp})
            nss = mean(NormalizedSquaredSignal_HJ(...
                [filtered_lfp.time{iRp}',filtered_lfp.lfp{iRp}]),2);
            time = filtered_lfp.time{iRp};
        else
            nss = cell(length(filtered_lfp.lfp{iRp}),1);
            for iW = 1:length(filtered_lfp.lfp{iRp})
                nss{iW} = mean(NormalizedSquaredSignal_HJ(...
                    [filtered_lfp.time{iRp}{iW}',filtered_lfp.lfp{iRp}{iW}]),2);
            end
            nss = cell2mat(nss);
            time = cell2mat(filtered_lfp.time{iRp}');
        end
        for jRp = 1:3 % the probe where ripples are detected from            
            if jRp<3
                rippletime = mat2cell(CA1_ripple_classified.(refList{jRp})(:,[1,3]),...
                    ones(size(CA1_ripple_classified.(refList{jRp}),1),1),2);
            else
                rippletime = mat2cell(CA1_ripple_classified.(refList{iRp+2})(:,[1,3]),...
                    ones(size(CA1_ripple_classified.(refList{iRp+2}),1),1),2);
            end
            inripple = sum(cell2mat(cellfun(@(x) time>=x(1) & time<=x(2),...
                rippletime,'UniformOutput',false)));
            if isempty(bincount{iRp,jRp})
                [bincount{iRp,jRp},bincount{iRp,jRp}] = deal(NaN(nS,length(binrange),2));
            end
            
            bincounttemp = hist(nss(inripple==0),binrange);
            bincounttemp(end) = sum(nss(inripple==0)>=binrange(end));
            bincount{iRp,jRp}(iS,:,1) = hist(nss(inripple==0),binrange)/sum(inripple==0); %in ripple
            bincount{iRp,jRp}(iS,:,2) = hist(nss(inripple>0),binrange)/sum(inripple>0); % in no ripple
        end
    end
end

%%
close all
stl = {':','-'};
titleList = {'dCA1 ripple';'iCA1 ripple';'Global ripple'};

ct = cbrewer('div','RdBu',7);
ct(2:6,:) = [];

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.45 20*0.14]);
for jRp = 1:3
   axes('Position',axpt(3,1,jRp,1,axpt(1,10,1,1:9)));
    hold on;
    for iRp = 1:2
        m = squeeze(mean(cumsum(bincount{iRp,jRp},2),1))';
        s = squeeze(std(cumsum(bincount{iRp,jRp},2),[],1))'/sqrt(nS);
        for i = 1:2
            when99 = find(m(i,:)>0.99,1,'first');
            fill([binrange(1:when99) flip(binrange(1:when99))],[m(i,1:when99)+s(i,1:when99),...
                flip(m(i,1:when99)-s(i,1:when99))],ct(iRp,:),'EdgeColor','none');
            plot(binrange(1:when99),m(i,1:when99),'Color',ct(iRp,:),'LineStyle',stl{i});
        end
    end
    alpha(0.2);
    set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'XTick',0:10:50,...
        'YTick',0:0.2:1,'XLim',[-1 50],'YLim',[0 1.01]);
    if jRp>1
        set(gca,'YTickLabel',[]);
    else
        ylabel('Cumulative fraction');
        xlabel('Normalized ripple band power');
    end
    
    title(titleList{jRp});
end
cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig1');
print(fHandle,'-depsc','-painters','ripple_power_cumsum_until99.ai');