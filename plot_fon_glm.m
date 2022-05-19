clearvars; clc; close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

target = [tag.area.hippo&tag.celltype.rs tag.area.hippo&tag.celltype.fs,...
    tag.area.vis&tag.celltype.rs tag.area.vis&tag.celltype.fs,...
    tag.area.thalamus tag.area.midbrain];
nTg = size(target,2);

refList = {'mrp';'lrp'};
[pval,beta] = deal(cell(nTg,2));
[area,ccf] = deal(cell(nTg,1));

for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_ripple_modulated.mat'],'fr_ripple_glm');
    
    for iTg = 1:nTg
        in = ismember(fr_ripple_glm.unit_id,tag.info.unit_id(target(:,iTg)));
        [~,idx] = ismember(fr_ripple_glm.unit_id(in),tag.info.unit_id);
        area{iTg} = [area{iTg};tag.info.structure(idx)];
        ccf{iTg} = [ccf{iTg};tag.info.ccf(idx)];
        for iRef = 1:2
            pval{iTg,iRef} = [pval{iTg,iRef}; fr_ripple_glm.(['p_',refList{iRef}])(in,:)];
            beta{iTg,iRef} = [beta{iTg,iRef}; fr_ripple_glm.(['b_',refList{iRef}])(in,:)];
        end
    end
end

time = -2:0.05:2;
clr = {'r';'b'};
tgList = {'Hippocampus-RS';'Hippocampus-FS';'Visual cortex-RS';'Visual cortex-FS';...
    'Thalamus';'Midbrain'};
refList = {'Dorsal';'Intermediate'};
x = {[1:3],[4:6],[8:10],[11:13],[15:17],[19:21]};
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 13 6]);
for iTg = 1:nTg
    axes('Position',axpt(21,11,x{iTg},1:3,[],[0.02 0.05]));
    hold on;
        plot([0 0],[0 1],'Color',[0.8 0.8 0.8]);
    for iRef = 1:2
        plot(time,nanmean(pval{iTg,iRef}<0.05),'Color',clr{iRef});
    end
    set(gca,'XTick',-2:2:2,'XTickLabel',[],'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35);
    if iTg<3
        ylim([0 0.9]);
        text(-1.7,0.8*0.97,['n = ',num2str(size(pval{iTg,iRef},1))],'FontSize',4);
        set(gca,'YTick',0:0.4:0.8,'YTickLabel',0:40:80);
    else
        ylim([0 0.45])
        text(-1.7,0.4*0.97,['n = ',num2str(size(pval{iTg,iRef},1))],'FontSize',4);
        set(gca,'YTick',0:0.2:0.4,'YTickLabel',0:20:40);
    end
    if ismember(iTg,[2 4])
        set(gca,'YTickLabel',[]);
    end
    if iTg==1
        ylabel('FON (%)','FontSize',5);
    end
    title(tgList{iTg},'FontSize',5);
    xlim([-2 2]);
    
    
    for iRef = 1:2
        axes('Position',axpt(21,11,x{iTg},[5:7]+(iRef-1)*3,[],[0.02 0.05]));
        hold on;
                plot([0 0],[0 1],'Color',[0.8 0.8 0.8]);

        plot(time,nanmean(pval{iTg,iRef}<0.05 & beta{iTg,iRef}>0),'Color',clr{iRef},'LineStyle','-');
        plot(time,nanmean(pval{iTg,iRef}<0.05 & beta{iTg,iRef}<0),'Color',clr{iRef},'LineStyle',':');
        set(gca,'XTick',-2:2:2,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35);
        if iTg<3
            ylim([0 0.8]);
            set(gca,'YTick',0:0.4:0.8,'YTickLabel',0:40:80);
        else
            ylim([0 0.4])
            set(gca,'YTick',0:0.2:0.4,'YTickLabel',0:20:40);
        end
        if ismember(iTg,[2 4])
            set(gca,'YTickLabel',[]);
        end
        if iTg==1
            ylabel({refList{iRef};'FON (%)'},'FontSize',5);
        end
        if iRef==1
            set(gca,'XTickLabel',[]);
        else
            xlabel('Time (s)','FontSize',5);
        end
        xlim([-2 2]);
    end
end
print(fHandle,'-dtiff','-r600','D:\OneDrive - University of California, San Francisco\figures\allen\fon\fon_glm.tif');

%% 
clr = {'b';'r'};
iTg = 5;
tgThal = {'LGd';'LP';'MG';'VPM'};
nThal = length(tgThal);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 6]);
for iThal = 1:nThal
    axes('Position',axpt(4,11,iThal,1:3,[],[0.02 0.05]));
    hold on;
        plot([0 0],[0 1],'Color',[0.8 0.8 0.8]);
        in = contains(area{iTg},tgThal{iThal});
    for iRef = 1:2
        plot(time,nanmean(pval{iTg,iRef}(in,:)<0.05),'Color',clr{iRef});
    end
    set(gca,'XTick',-2:2:2,'XTickLabel',[],'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35);
    if iThal<3
        ylim([0 0.6]);
        text(-1.7,0.6*0.9,['n = ',num2str(sum(in))],'FontSize',4);
        set(gca,'YTick',0:0.3:0.6,'YTickLabel',0:30:60);
    else
        ylim([0 0.6])
        text(-1.7,0.6*0.9,['n = ',num2str(sum(in))],'FontSize',4);
        set(gca,'YTick',0:0.3:0.6,'YTickLabel',0:30:60);
    end
    if iThal==1
        ylabel('FON (%)','FontSize',5);
    else
        set(gca,'YTickLabel',[]);
    end
    title(tgThal{iThal},'FontSize',5);
    xlim([-2 2]);
    
    
    for iRef = 1:2
        axes('Position',axpt(4,11,iThal,[5:7]+(iRef-1)*3,[],[0.02 0.05]));
        hold on;
                plot([0 0],[0 1],'Color',[0.8 0.8 0.8]);

        plot(time,nanmean(pval{iTg,iRef}(in,:)<0.05 & beta{iTg,iRef}(in,:)>0),'Color',clr{iRef},'LineStyle','-');
        plot(time,nanmean(pval{iTg,iRef}(in,:)<0.05 & beta{iTg,iRef}(in,:)<0),'Color',clr{iRef},'LineStyle',':');
        set(gca,'XTick',-2:2:2,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35);

        if iThal==1
            ylabel({refList{iRef};'FON (%)'},'FontSize',5);
            set(gca,'YTick',0:0.3:0.6,'YTickLabel',0:30:60);
        else
            set(gca,'YTickLabel',[]);
        end
        if iRef==1
            set(gca,'XTickLabel',[]);
        else
            xlabel('Time (s)','FontSize',5);
        end
        xlim([-2 2]);
        ylim([0 0.6]);
    end
end
print(fHandle,'-dtiff','-r600','D:\OneDrive - University of California, San Francisco\figures\allen\fon\fon_glm_thalamus.tif');
