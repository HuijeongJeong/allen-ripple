clearvars; clc; close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nSession = length(sessionList);

fName = {'ripple'};
fRange = [100 250];

structureList = {'VISp';'VISl';'VISal';'VISrl';'VISpm';'VISam'};
nStructure = length(structureList);
totalNss = cell(nStructure,3);

for iS = 1:nSession
    iS
    sessionDir = ['D:\OneDrive\1.allen-andermann\session_',num2str(sessionList(iS))];
    try
        power_ripple =  load([sessionDir,'\session_',num2str(sessionList(iS)),...
            '_ripple_modulated_classified.mat'],'ripple_power_ripple');
    catch
        continue;
    end
    structure = power_ripple.ripple_power_ripple.structure;
    nss = cellfun(@nanmean,power_ripple.ripple_power_ripple.nss,'UniformOutput',false);
    time = power_ripple.ripple_power_ripple.time;
    nss = mat2cell(cell2mat(nss),size(nss,1),repmat(length(time),1,3));
    for iStr = 1:nStructure
        inStr = strcmp(structure,structureList{iStr});
        if sum(inStr)==0
            continue; end
        totalNss(iStr,:) = cellfun(@(y,z) [y;z],squeeze(totalNss(iStr,:)),...
            cellfun(@(x) nanmean(x(inStr,:)),nss,'UniformOutput',false),'UniformOutput',false);
    end
end

%%
dir = 'D:\OneDrive - University of California, San Francisco\figures\allen\figS4';

%% plot average nss
clr = {'k';'r';'b';'k'};
clr_dim = {[0.6 0.6 0.6],[1 0.6 0.6],[0.6 0.6 1],[0.6 0.6 0.6]};
ylabelList = {{'Global ripple';'Norm. ripple power'};{'dCA1 ripple';'Norm. ripple power'};...
    {'iCA1 ripple';'Norm. ripple power'};{'\DeltaNorm. ripple power';'(dCA1-iCA1)'}};
structure_acronym = {'V1';'LM';'AL';'RL';'PM';'AM'};

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 15 8]);
for iRef = 1:4
    for iArea = 1:nStructure
        if iRef<4
            axes('Position',axpt(nStructure,17,iArea,[1:4]+(iRef-1)*4,[],[0.025 0.05]))
        data = totalNss{iArea,iRef};
        plot([0 0],[-1 6],'k:','LineWidth',0.35);
        else
            axes('Position',axpt(nStructure,17,iArea,14:17,[],[0.025 0.05]))
            data = totalNss{iArea,2}-totalNss{iArea,3};
            
            plot([0 0],[-4 4],'k:','LineWidth',0.35);
        end
        hold on;
        m = nanmean(data);
%         s = nanstd(data)/sqrt(size(data,1));
%         fill([time flip(time)],[m+s flip(m-s)],clr{iRef},'EdgeColor','none');
        plot(time,data,'Color',clr_dim{iRef},'LineWidth',0.35);
        plot(time,m,'Color',clr{iRef},'LineWidth',1);
        alpha(0.2)
        xlim([-0.2 0.2]);
        set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
            'XTick',-0.2:0.2:0.2,'YTick',-3:3:6);
        if iRef<4
            ylim([-0.5 6]);
            if iRef<3
                set(gca,'XTickLabel',[]);
                if iRef==1
                   title([structure_acronym{iArea},' (n=',...
                       num2str(size(data,1)),')'],'FontSize',8);
                end
            end
        else
            ylim([-4 4])
            if iArea==1
            xlabel('Time from ripple onset (s)','FontSize',8);
            end
        end
        if iArea==1
            ylabel(ylabelList{iRef},'FontSize',8);
        else
            set(gca,'YTickLabel',[]);
        end
    end
end
print(fHandle,'-depsc','-painters',[dir,'\average_nss.ai']);

%% plot area under curve
close all
refName = {'Global';'dCA1';'iCA1'};
ylimit = [3 1.75 1.75];
bin = 0.1;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 12 4]);
% auc_pre = cellfun(@(x) trapz(time(time>=-bin&time<0),x(:,time>=-bin& time<0)')/bin,...
%     totalNss,'UniformOutput',false);
% auc_post = cellfun(@(x) trapz(time(time>0&time<=bin),x(:,time>0 & time<=bin)')/bin,...
%     totalNss,'UniformOutput',false);

auc_pre = cellfun(@(x) sum(x(:,time>=-bin& time<0),2)*mean(diff(time)),...
    totalNss,'UniformOutput',false);
auc_post = cellfun(@(x) sum(x(:,time>0 & time<=bin),2)*mean(diff(time)),...
    totalNss,'UniformOutput',false);
for iRef = 1:4
    if iRef<4
    axes('Position',axpt(13,1,[1:3]+(iRef-1)*3,1,axpt(10,10,2:10,1:9)));
    else
    axes('Position',axpt(13,1,11:13,1,axpt(10,10,2:10,1:9)));
    end
        hold on;
        plot([0 nStructure+1],[0 0],'k:','LineWidth',0.35);
    for i = 1:4        
        if iRef<4
            m_pre = cellfun(@nanmean,auc_pre(:,iRef));
            s_pre = cellfun(@(x) nanstd(x)/sqrt(length(x)),auc_pre(:,iRef));
            m_post = cellfun(@nanmean,auc_post(:,iRef));
            s_post = cellfun(@(x) nanstd(x)/sqrt(length(x)),auc_post(:,iRef));
        else
            m_pre = cellfun(@(x,y) nanmean(x-y),auc_pre(:,2),auc_pre(:,3));
            s_pre = cellfun(@(x,y) nanstd(x-y)/sqrt(length(x)),auc_pre(:,2),auc_pre(:,3));
            m_post = cellfun(@(x,y) nanmean(x-y),auc_post(:,2),auc_post(:,3));
            s_post = cellfun(@(x,y) nanstd(x-y)/sqrt(length(x)),auc_post(:,2),auc_post(:,3));
        end
        
        errorbar(1:nStructure,m_pre,s_pre,'Color','k','LineWidth',0.5,...
            'LineStyle',':','CapSize',4);
        errorbar(1:nStructure,m_post,s_post,'Color','k','LineWidth',0.5,'CapSize',4);
        xlim([0 nStructure+1]);
    end
    set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
            'XTick',1:nStructure,'XTickLabel',structure_acronym,'XTickLabelRotation',45)
    if iRef<4
        ylim([-0.01 0.1]);
        set(gca,'YTick',0:0.05:0.1);
        if iRef>1
            set(gca,'YTickLabel',[]);
        else
            ylabel('AUC');
        end
    else
        [h,p,~,stat] = cellfun(@(x,y) ttest(x,y),cellfun(@(x,y) x-y,auc_pre(:,2),auc_pre(:,3),'UniformOutput',false),...
            cellfun(@(x,y) x-y,auc_post(:,2),auc_post(:,3),'UniformOutput',false),'UniformOutput',false);
        for iS = 1:nStructure
            if h(iS)==1
               text(iS,0.05,'*'); 
            end
        end
        ylim([-0.05 0.05]);
        set(gca,'YTick',-0.05:0.05:0.05);
        ylabel({'\DeltaAUC';'(dCA1-iCA1)'},'FontSize',8);
        
    end
end
print(fHandle,'-depsc','-painters',[dir,'\ripple_power_auc.ai']);



