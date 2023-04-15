clearvars;clc;close all;
rng(5)

load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');

%% setting up 
k = 3; % number of cluster 
nIter = 100;

invis = tag.area.vis6;
session = tag.info.session_id(invis);
sessionlist = unique(session);
celltype = [tag.celltype.rs(invis),tag.celltype.fs(invis)];
area = tag.info.structure(invis);
[in,idx] = ismember(tag.info.unit_id(invis),unit_id.vis);
cidx = zeros(sum(invis),1);
cidx(in) = cluster_idx.vis{k-1}(idx(in));

areaList = {'VISp';'VISl';'VISal';'VISrl';'VISam';'VISpm'};
nArea = length(areaList);
fon = cell(k,1);
iCT = 1;
for iA = 1:nArea
    cidx_area = cidx(strcmp(area,areaList{iA}) & celltype(:,iCT));
    cidx_area_sig = cidx(strcmp(area,areaList{iA}) & celltype(:,iCT) & cidx>0);
    
    for iIter = 1:nIter
        cidx_area_sample = randsample(cidx_area,100);
        cidx_area_sig_sample = randsample(cidx_area_sig,100);
        for iClst = 1:k
            if isempty(fon{iClst,1})
                [fon{iClst,1},fon{iClst,2}] = deal(NaN(nIter,nArea));
            end
            fon{iClst,1}(iIter,iA) = sum(cidx_area_sample==iClst);
            fon{iClst,2}(iIter,iA) = sum(cidx_area_sig_sample==iClst);
        end
    end
end


%%
x = [repmat(1:nArea,2,1); NaN(1,nArea)];

ct = cbrewer('qual','Dark2',3);

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20/3 20*0.15]);
xx = repmat([2:nArea],100,1);
for iClst = 1:k
h(iClst) = axes('Position',axpt(3,1,iClst,1)); hold on;
fon_sig = sort(fon{iClst,2});
y = [fon_sig(5,:); fon_sig(95,:); NaN(1,nArea)];
plot(x(:),y(:),'Color',ct(iClst,:));
scatter(x(1,:),mean(fon_sig),2,ct(iClst,:),'s','filled')
title(['Cluster ',num2str(iClst)],'Color',ct(iClst,:));
[~,p(iClst),~,stat] = ttest2(sum(fon{iClst,2}(:,2:3),2),...
    sum(fon{iClst,2}(:,end-1:end),2));
t(iClst) = stat.tstat;
if iClst==1
    ylabel('FON (%)');
end
end
set(h,'XLim',[0.5 6.5],'YLim',[0 65],'Box','off','TickDir','out','FontSize',7,...
    'YTick',0:20:60,'XTick',1:nArea,'XTickLabel',...
    {'V1';'LM';'AL';'RL';'AM';'PM'},'XTickLabelRotation',45)
set(h(2:end),'YTickLabel',[])
% cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig1')
% print(fHandle,'-depsc','-painters','fon.ai')

%% 
fon_s = nan(length(sessionlist),2);
for iS = 1:length(sessionlist)
    if sum(cidx>0&celltype(:,iCT)&session==sessionlist(iS))==0
        continue;
    end
    ins = celltype(:,iCT)&session==sessionlist(iS);
%     for iclst = 1:4
%         fon_s(iS,iclst,1) = sum(cidx==iclst-1&ins&ismember(area,{'VISl','VISal'}));
%         fon_        
%     end
    fon_s(iS,1) = sum(cidx>0&ins&ismember(area,{'VISl','VISal'}))/sum(ins&ismember(area,{'VISl','VISal'}));
    fon_s(iS,2) = sum(cidx>0&ins&ismember(area,{'VISam','VISpm'}))/sum(ins&ismember(area,{'VISam','VISpm'}));
end
fon_s(sum(isnan(fon_s),2)==2,:) = [];
% fon_s = rmmissing(fon_s);
[~,p,~,stat] = ttest(fon_s(:,1),fon_s(:,2));

%%

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.5 20*0.17]);
hold on;
fon_sig = fon{1,1}+fon{2,1}+fon{3,1};
[~,pp] = ttest2(sum(fon_sig(:,2:3),2),sum(fon_sig(:,end-1:end),2));
fon_sig = sort(fon_sig);
y = [fon_sig(5,:); fon_sig(95,:); NaN(1,nArea)];
plot(x(:),y(:),'k');
scatter(x(1,:),mean(fon_sig),2,'ks','filled')
ylabel('% ripple modulated')
set(gca,'XLim',[0.5 6.5],'YLim',[0 40],'Box','off','TickDir','out','FontSize',7,...
    'YTick',0:20:40,'XTick',1:nArea,'XTickLabel',...
    {'V1';'LM';'AL';'RL';'AM';'PM'},'XTickLabelRotation',45)
cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig1')
print(fHandle,'-depsc','-painters','fon_ripplemodulated.ai')
