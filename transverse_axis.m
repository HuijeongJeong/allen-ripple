clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

distance2medial = cell(nS,6);
coordinates = cell(nS,6);
for iS = 1:nS
    iS
%     load([sdir(sessionList(iS)),'_cellTable.mat'],'T');
    load([sdir(sessionList(iS)),'_ripples.mat'],'total_CA1_ripple');
    load([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_CA1_ripple');
    
    n = size(spontaneous_CA1_ripple,1);
    medialidx = find(strcmp(spontaneous_CA1_ripple.relative_location,'medial'));
    lateralidx = find(strcmp(spontaneous_CA1_ripple.relative_location,'lateral'));
    ccf = cell2mat(spontaneous_CA1_ripple.ccf_coordinate);
    d2m = abs(ccf-ccf(medialidx,:));
    [~,sortidx] = sort(sum(d2m,2));
    distance2medial(iS,1:n) = mat2cell(d2m(sortidx,:),ones(n,1),3);
    
    coordinates(iS,1:n) = total_CA1_ripple.ccf_coordinate(sortidx);
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
k = boundary(data(:,1),data(:,2),data(:,3),0.95);

%%
clr = [0,0,0;cbrewer('div','RdBu',9)];
clr(4:7,:) = [];
nprobes = sum(cellfun(@(x) ~isempty(x),coordinates),2);
inanal = nprobes>=5;
%%
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 4]);
for i = 1:2
    axes('Position',axpt(2,1,i,1,axpt(1,10,1,1:9),[0.18 0.05]));
    hold on;
    trisurf(k,data(:,1),data(:,2),data(:,3),...
        'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.5,'EdgeColor','none');
    %%
    for irp = [1,2,3,5]
        if irp<4
            data_sub = coordinates(inanal,irp);
%             dis_sub = distance2medial(inanal,irp);
        else
            data_sub = cellfun(@(x,y) squeeze(x(:,y-(5-irp),:))',...
                mat2cell(coordinates(inanal,:),ones(sum(inanal),1),6),...
                num2cell(nprobes(inanal)),'UniformOutput',false);
            data_sub = cat(1,data_sub{:});
%             
%             dis_sub = cellfun(@(x,y) squeeze(x(:,y-(5-irp),:))',...
%                 mat2cell(distance2medial(inanal,:),ones(sum(inanal),1),6),...
%                 num2cell(nprobes(inanal)),'UniformOutput',false);
%             dis_sub = cat(1,dis_sub{:});
        end
       
        data_sub = cell2mat(data_sub)/1000;
        scatter3(data_sub(:,1),data_sub(:,3),data_sub(:,2),2,clr(irp,:),'filled');
        scatter3(mean(data_sub(:,1)),mean(data_sub(:,3)),mean(data_sub(:,2)),...
            8,clr(irp,:),'filled');
%         if i==1
%         mean(data_sub)
%         mean(dis_sub)
%         end
    end
    %%
    
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
print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\figures\2.allen\revision\FigS10\location_of_probes_2.ai');
%%

