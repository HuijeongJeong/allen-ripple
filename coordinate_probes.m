clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

coordinates = cell(nS,1);
for iS = 1:nS
    iS
%     load([sdir(sessionList(iS)),'_cellTable.mat'],'T');
    load([sdir(sessionList(iS)),'_ripples.mat'],'total_CA1_ripple');
    
    coordinates{iS} = sortrows(cell2mat(total_CA1_ripple.ccf_coordinate),3);
end

%%
close all
meshDir = 'D:\OneDrive\1.allen-andermann\mouse_ccf\annotation\ccf_2017\structure_meshes\';
url = 'http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/structure_meshes/';
fileID = fopen([meshDir,'CA1_382.txt']);
if fileID==-1
    websave([meshDir,'CA1_382.txt'],[url,'/382.obj']);
    fileID = fopen([meshDir,'CA1_382.txt']);
end
C = textscan(fileID,'%s %f %f %f','HeaderLines',1);
in = strcmp(C{1}(1:end-1),'v');

data = cell2mat(cellfun(@(x) x(in),C([2 4 3]),'UniformOutput',false))/1000;

medial = cell2mat(cellfun(@(x) x(1,:)/1000,coordinates,'UniformOutput',false));
lateral = cell2mat(cellfun(@(x) x(end,:)/1000,coordinates,'UniformOutput',false));
other = cell2mat(cellfun(@(x) x(2:end-1,:)/1000,coordinates,'UniformOutput',false));

k = boundary(data(:,1),data(:,2),data(:,3),0.95);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 4]);
for i = 1:2
    axes('Position',axpt(2,1,i,1,axpt(1,10,1,1:9),[0.18 0.05]));
    hold on;
    trisurf(k,data(:,1),data(:,2),data(:,3),...
        'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.5,'EdgeColor','none');
    
    scatter3(other(:,1),other(:,3),other(:,2),5,[0.3 0.3 0.3],'filled');
    scatter3(medial(:,1),medial(:,3),medial(:,2),5,'r','filled');
    scatter3(lateral(:,1),lateral(:,3),lateral(:,2),5,'b','filled');
    if i==1
        view([90 0]);
    elseif i==2
        view([90 90]);
    end
    set(gca,'ZDir','reverse','YLim',[5.5 10],'XLim',[5.5 10],'ZLim',[1 7],'Box','off','TickDir',...
        'out','FontSize',7,'LineWidth',0.35,'XTick',6:2:8,'YTick',6:2:10,'ZTick',2:2:6);
    xlabel('Anterior-posterior (mm)')
    ylabel('Medial-lateral (mm)')
    zlabel('Dorsal-ventral (mm)')
end
%%
print(fHandle,'-dtiff','-r600','D:\OneDrive - UCSF\data\allen\revision\location_of_all_probes.tif');

%%
close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5 4]);
data = lateral(:,[1,3,2])-medial(:,[1,3,2]);
hold on;
for i = 1:3
    bar(i,mean(data(:,i)),0.8,'FaceColor',[0.8 0.8 0.8]);
    scatter(rand([1,nS])*0.6-0.3+i,data(:,i),5,[0.4 0.4 0.4],'filled');
    errorbar(i,mean(data(:,i)),std(data(:,i))/sqrt(nS),'k');
end

data = sqrt(sum((lateral-medial).^2,2));
bar(4.5,mean(data),0.8,'FaceColor',[0.8 0.8 0.8]);
scatter(rand([1,nS])*0.6-0.3+4.5,data,5,[0.4 0.4 0.4],'filled');
errorbar(4.5,mean(data),std(data)/sqrt(nS),'k');
ylim([0 3]);
xlim([0 5.5]);
set(gca,'Box','off','Tickdir','out','FontSize',8,'XTick',[1,2,3,4.5],...
    'XTickLabel',{'AP';'ML';'DV';'total'},'XTickLabelRotation',45,'YTick',0:1:3);
ylabel('\Delta distance (mm)');
print(fHandle,'-dtiff','-r600','D:\OneDrive - UCSF\data\allen\revision\distance_bw_probes.tif');
