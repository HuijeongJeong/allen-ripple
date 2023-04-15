clearvars; clc; close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
annotation = readtable('D:\OneDrive\1.allen-andermann\annotation.csv');
sessionList = session_metric.session_id(session_metric.immobile_period>100 & session_metric.distance_bw_ml_probes>2000);
nSession = length(sessionList);


%% load CA1 structure mesh from CCFv 
meshDir = 'D:\OneDrive\1.allen-andermann\mouse_ccf\annotation\ccf_2017\structure_meshes\';
url = 'http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/structure_meshes/';

areaList = {'VISp';'VISl';'VISrl';'VISal';'VISpm';'VISam'};
areaid = cellfun(@(x) annotation.id(strcmp(x,annotation.acronym)),areaList);
nA = length(areaList);
%%

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20/8*2 20/8*2]);
hold on;

for iA = 1:nA
fileID = fopen([meshDir,areaList{iA},'_',num2str(areaid(iA)),'.txt']);
if fileID==-1
    websave([meshDir,areaList{iA},'_',num2str(areaid(iA)),'.txt'],[url,'/',num2str(areaid(iA)),'.obj']);
    fileID = fopen([meshDir,areaList{iA},'_',num2str(areaid(iA)),'.txt']);
end
C = textscan(fileID,'%s %f %f %f','HeaderLines',1);
in = strcmp(C{1}(1:end-1),'v') & C{4}>5000;

data = cell2mat(cellfun(@(x) x(in),C([2 4 3]),'UniformOutput',false))/1000;

k = boundary(data(:,1),data(:,2),data(:,3),0.95);
trisurf(k,data(:,1),data(:,2),data(:,3),...
    'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.5,'EdgeColor','none');
end
ylim([6.5 10.5])
xlim([6.5 10.5]);
zlim([0 3]);
view([80 60])
set(gca,'ZDir','reverse','Box','off','TickDir','out','FontSize',5,...
    'LineWidth',0.35,'XTick',7:10,'YTick',7:10,'ZTick',0:3,...
    'XTickLabel',7:10,'YTickLabel',7:10,'ZTickLabel',0:3);

%%
for iS = 1:nSession
    sessionDir = ['D:\OneDrive\1.allen-andermann\session_',num2str(sessionList(iS))];
    try
        load([sessionDir,'\session_',num2str(sessionList(iS)),'_ripples.mat'],...
            'spontaneous_CA1_ripple');
    catch
        continue;
    end
    
    ccf = cell2mat(spontaneous_CA1_ripple.ccf_coordinate)/1000;
    loc_idx = spontaneous_CA1_ripple.relative_location;

    m = cellfun(@(x) strcmp(x,'medial'),loc_idx);
    l = cellfun(@(x) strcmp(x,'lateral'),loc_idx);
    scatter3(ccf(m,1),ccf(m,3),ccf(m,2),5,'r','filled'); % most medially located probe
    scatter3(ccf(l,1),ccf(l,3),ccf(l,2),5,'b','filled'); % most medially located probe
end
set(gca,'ZDir','reverse','YLim',[5.5 10],'Box','off','TickDir',...
    'out','FontSize',7,'LineWidth',0.35,'XTick',6:2:8,'YTick',6:2:10,'ZTick',2:2:6);
view([80 25])
xlabel('Anterior-posterior (mm)')
ylabel('Medial-lateral (mm)')
zlabel('Dorsal-ventral (mm)')
cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig1')
print('-depsc','-painters','probe_loc.ai');
