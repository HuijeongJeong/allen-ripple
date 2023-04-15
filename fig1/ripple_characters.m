clearvars; clc; close all

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

refList = {'medial';'lateral';'global'};
globalList = {'global_m';'global_l'};
typeList = {'rs';'fs'};
nclst = 3;

[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
session_id = tag.info.session_id(idx);

unitid = unit_id.vis;
sessionList = unique(session_id);
nS = length(sessionList);
binsize = 0.05;
win = [-0.5 0.5];

nIter = 1000;
ripplecorr = NaN(nS,diff(win)/binsize);
ripplecorr_sf = NaN(nS,nIter,diff(win)/binsize);

[ripplelength,ripplestrength,gripplefraction] = deal(NaN(nS,5));
gripple_mstartfraction = nan(nS,1);
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_Events.mat'],'probe');
    load([sdir(sessionList(iS)),'_ripples.mat'],'CA1_ripple_classified','spontaneous_win','spontaneous_anal_win');
    
    ripplehist = cell2mat(cellfun(@(x) histc(CA1_ripple_classified.(refList{2})(:,1),...
        [win(1):binsize:win(2)]+x)',num2cell(CA1_ripple_classified.(refList{1})(:,1)),...
        'UniformOutput',false));
    ripplecorr(iS,:) = mean(ripplehist(:,1:end-1));
    for iIter = 1:nIter
        shuffled = cell(size(spontaneous_anal_win,1),1);
        for iw = 1:size(spontaneous_anal_win,1)
            inwindow = CA1_ripple_classified.(refList{2})(:,1)>=spontaneous_anal_win(iw,1) &...
                CA1_ripple_classified.(refList{2})(:,1)<=spontaneous_anal_win(iw,2);
            jitter = rand([sum(inwindow),1])*diff(spontaneous_anal_win(iw,:))-...
                0.5*diff(spontaneous_anal_win(iw,:));
            
            shuffled_temp = CA1_ripple_classified.(refList{2})(inwindow,1)+jitter-spontaneous_anal_win(iw,1);
            shuffled_temp(shuffled_temp<0) = shuffled_temp(shuffled_temp<0)+diff(spontaneous_anal_win(iw,:));
            shuffled_temp(shuffled_temp>diff(spontaneous_anal_win(iw,:))) =...
                shuffled_temp(shuffled_temp>diff(spontaneous_anal_win(iw,:)))-diff(spontaneous_anal_win(iw,:));
            
            shuffled{iw} = sort(shuffled_temp+spontaneous_anal_win(iw,:));
        end
        ripplehist_sf = cell2mat(cellfun(@(x) histc(cell2mat(shuffled),...
            [win(1):binsize:win(2)]+x)',num2cell(CA1_ripple_classified.(refList{1})(:,1)),...
            'UniformOutput',false));
        ripplecorr_sf(iS,iIter,:) = mean(ripplehist_sf(:,1:end-1));
    end
    for iRp = 1:3
        if iRp<3
            ripplelength(iS,iRp) = mean(diff(CA1_ripple_classified.(refList{iRp})(:,[1 3]),[],2));
            ripplestrength(iS,iRp) = mean(CA1_ripple_classified.(refList{iRp})(:,4));
            gripplefraction(iS,iRp) = mean(CA1_ripple_classified.([refList{iRp},'_overlapInd']));
        else
            ripplelength(iS,iRp) = mean(diff(CA1_ripple_classified.(refList{iRp}),[],2));
            ripplestrength(iS,iRp) = mean(max(CA1_ripple_classified.global_strength,[],2));
            gripple_startm_ind = ismember(CA1_ripple_classified.(refList{iRp})(:,1),...
                CA1_ripple_classified.global_m(:,1));            
            for i = 1:2
               ripplelength(iS,iRp+i) = mean(diff(CA1_ripple_classified.(refList{iRp})(gripple_startm_ind==2-i,:),...
                   [],2));
               ripplestrength(iS,iRp+i) = mean(max(CA1_ripple_classified.global_strength(gripple_startm_ind==2-i,:),[],2));
               gripple_mstartfraction(iS) = mean(gripple_startm_ind);
            end
        end
        
    end
end

% %%
% fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2 20*0.17]);
% hold on;
% plot(1:2,ripplelength*1000,'Color',[0.6 0.6 0.6],'linewidth',0.35);
% errorbar(1:2,mean(ripplelength*1000),std(ripplelength*1000)/sqrt(nS),'k','CapSize',3,'LineWidth',1);
% xlim([0.5 2.5])
% ylim([20 60]);
% set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'TickDir','out',...
%     'XTick',1:2,'XTickLabel',{'dCA1';'iCA1'},'XTickLabelRotation',45,'YTick',20:20:60);
% ylabel('Ripple length (ms)');
% print(fHandle,'-depsc','-painters','ripple_length.ai');

%%
dir = 'D:\OneDrive - University of California, San Francisco\figures\allen\figS1';
%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5 4]);
h(1) = axes('Position',axpt(2,1,1,1,axpt(10,10,2:10,1:8),[0.2 0.05])); 
hold on;
plot(1:3,ripplelength(:,1:3)*1000,'Color',[0.6 0.6 0.6]);
errorbar(1:3,mean(ripplelength(:,1:3)*1000),std(ripplelength(:,1:3)*1000)/sqrt(nS),'k','CapSize',3);
% [stats,rm] = simple_mixed_anova(ripplelength(:,1:3));
% multcompare(rm,'WS01');

h(2) = axes('Position',axpt(2,1,2,1,axpt(10,10,2:10,1:8),[0.2 0.05])); 
hold on;
plot(1:3,ripplestrength(:,1:3),'Color',[0.6 0.6 0.6]);
errorbar(1:3,mean(ripplestrength(:,1:3)),std(ripplestrength(:,1:3))/sqrt(nS),'k','CapSize',3);
set(h,'XLim',[0.5 3.5],'XTick',1:3,'XTickLabel',{'dCA1';'iCA1';'Global'},...
    'XTickLabelRotation',45,'Box','off','TickDir','out','FontSize',8);
set(h(2),'YLim',[0 100],'YTick',0:50:100);
set(h(1),'YLim',[20 100],'YTick',20:40:100);
ylabel(h(1),'Ripple length (ms)');
ylabel(h(2),'Ripple strength (NSS)');
% [~,rm] = simple_mixed_anova(ripplestrength(:,1:3));
% multcompare(rm,'WS01');
% print(fHandle,'-depsc','-painters',[dir,'\dca1_ica1_ripple_comparison.ai']);

%%
p = [];
t = [];
for iData = 1:2
    for i = 1:2
        for j = i+1:3
            switch iData
                case 1
                    [~,ptemp,~,stat] = ttest(ripplelength(:,i),ripplelength(:,j));
                case 2
                    [~,ptemp,~,stat] = ttest(ripplestrength(:,i),ripplestrength(:,j));
            end
            p = [p;ptemp*6];
            t = [t; stat.tstat];
        end
    end
end

%%
close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 9 4]);
h(1) = axes('Position',axpt(11,1,1:3,1,axpt(10,10,2:10,1:8),[0.2 0.05])); 
hold on;
bar(0.5,mean(gripple_mstartfraction*100),1,'FaceColor',[0.6 0.6 0.6]);
errorbar(0.5,mean(gripple_mstartfraction*100),std(gripple_mstartfraction*100)/sqrt(nS),'k');
scatter(rand(nS,1)*0.8+0.1,gripple_mstartfraction*100,2,'k','filled');
plot([-0.5 1.5],[50 50],'k:')
[~,p,~,stat] = ttest(gripple_mstartfraction,0.5)

h(2) = axes('Position',axpt(11,1,4:7,1,axpt(10,10,2:10,1:8),[0.2 0.05])); 
hold on;
plot(1:2,ripplelength(:,4:5)*1000,'Color',[0.6 0.6 0.6]);
errorbar(1:2,mean(ripplelength(:,4:5)*1000),std(ripplelength(:,4:5)*1000)/sqrt(nS),'k','CapSize',3);
[~,p] = ttest(ripplelength(:,4),ripplelength(:,5))

h(3) = axes('Position',axpt(11,1,8:11,1,axpt(10,10,2:10,1:8),[0.2 0.05])); 
hold on;
plot(1:2,ripplestrength(:,4:5),'Color',[0.6 0.6 0.6]);
errorbar(1:2,mean(ripplestrength(:,4:5)),std(ripplestrength(:,4:5))/sqrt(nS),'k','CapSize',3);
set(h,'Box','off','TickDir','out','FontSize',8);
set(h(2:3),'XLim',[0.5 2.5],'XTick',1:2,'XTickLabel',{'dCA1 start';'iCA1 start'},...
    'XTickLabelRotation',45);
set(h(3),'YLim',[0 100],'YTick',0:50:100);
set(h(2),'YLim',[20 100],'YTick',20:40:100);
set(h(1),'YLim',[0 100],'YTick',0:50:100,'XLim',[-0.5 1.5],'XTick',[]);
ylabel(h(1),{'Fraction of global ripple'; 'starting from dCA1 (%)'});
ylabel(h(2),{'Global ripple';'length (ms)'});
ylabel(h(3),{'Global ripple'; 'strength (NSS)'});
[~,p] = ttest(ripplestrength(:,4),ripplestrength(:,5))
print(fHandle,'-depsc','-painters',[dir,'\global_ripple_detail.ai']);

%%
close all;
ripplecorr_sf_sorted = sort(ripplecorr_sf,2);

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 4 4]);
bar(win(1)+binsize/2:binsize:win(2),mean(ripplecorr),'FaceColor',[0.6 0.6 0.6],'EdgeColor','none');
hold on;
errorbar(win(1)+binsize/2:binsize:win(2),mean(ripplecorr),...
    std(ripplecorr)/sqrt(nS),'LineStyle','none','Color','k','CapSize',1)
plot(win(1)+binsize/2:binsize:win(2),mean(squeeze(ripplecorr_sf_sorted(:,0.95*1000,:))),'r--')
plot(win(1)+binsize/2:binsize:win(2),mean(squeeze(ripplecorr_sf_sorted(:,0.99*1000,:))),'r:')
plot([0 0],[0 0.05],'k:','LineWidth',0.35)
ylim([0 0.04])
xlim([-0.5 0.5])
ylabel({'Normalized';'cross-correlation'});
xlabel({'Time of iCA1 ripple';'from dCA1 ripple (s)'});
set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,'TickDir','out',...
    'XTick',-0.5:0.5:0.5,'YTick',0:0.02:0.04);
print(fHandle,'-depsc','-painters',[dir,'\ripple_cross_corr.ai']);

