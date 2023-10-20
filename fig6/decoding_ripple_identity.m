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

win = [-3 3];
bin = 0.5;
step = 0.1;
refList = {'medial','lateral'};
nIter = 100;

accuracy = cell(nS,5);

for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_win');
    load([sdir(sessionList(iS)),'_ripples.mat'],'CA1_ripple_classified');
    
    invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis & tag.celltype.rs));
    
    [in,idx] = ismember(T.unit_id,unitid);
    clstidx = zeros(length(T.unit_id),1);
    clstidx(in) = cluster_idx.vis{nclst-1}(idx(in));
    clstidx = clstidx(invis);
    clstidx = [clstidx==0, clstidx==1, clstidx==2, clstidx==3];
    
%     area{iS} = T.ecephys_structure_acronym(invis);

    spkhist = cell(sum(invis),2);
    for iRef = 1:2
        spkhist(:,iRef) = cellfun(@(y) cell2mat(cellfun(@(x) movsum(histc(y,[win(1):step:win(2)]+x)',...
            bin/step,'Endpoints','discard'),num2cell(CA1_ripple_classified.(refList{iRef})(:,1)),...
            'UniformOutput',false)),T.spike_time(invis),'UniformOutput',false);
    end
    
    nbin = size(spkhist{1,1},2); 
    nmin = min(sum(clstidx));
%     nmin = floor(min(sum(clstidx))/4);
    if nmin==0
        continue
    end
    ntraining = floor(min([size(CA1_ripple_classified.medial,1),...
        size(CA1_ripple_classified.lateral,1)])/2);
    trial = cellfun(@(x) [1:size(x,1)]',spkhist(1,:),'UniformOutput',false);
    
    
    for iIter = 1:nIter
        trainingtrial = cellfun(@(x) randsample(x,ntraining),trial,'UniformOutput',false);
        testtrial = cellfun(@(x,y) randsample(x(~ismember(x,y)),ntraining),trial,trainingtrial,'UniformOutput',false);
%         inclst = cell(1,4);
        for iclst = 1:4
            if isempty(accuracy{iS,iclst})
                accuracy{iS,iclst} = NaN(nIter,nbin);
            end
            inclst = randsample(find(clstidx(:,iclst)),nmin); 
%            inclst{1,iclst} = randsample(find(clstidx(:,iclst)),floor(nmin/4)*3); 
%            inclst{2,iclst} = randsample(find(clstidx(:,iclst)),nmin*4);
%         end
%         for i = 1:5
%             if i==1
%                 j = 1;
%             else
%                 j = 2;
%             end
                        
%             if isempty(accuracy{iS,i})
%                 accuracy{iS,i} = NaN(nIter,nbin);
%             end
            for iB = 1:nbin
                [trainingspk,testingspk] = deal(cell(1,2));
%                 [trainingspk,testingspk] = deal(cell(4,2));
                 for iRef = 1:2
                      trainingspk{iRef} = cell2mat(cellfun(@(x) x(trainingtrial{iRef},iB)',...
                        spkhist(inclst,iRef),'UniformOutput',false));
                    testingspk{iRef} = cell2mat(cellfun(@(x) x(testtrial{iRef},iB)',...
                        spkhist(inclst,iRef),'UniformOutput',false));
                    
%                     for iclst = 1:4
%                         
%                     trainingspk{iclst,iRef} = cell2mat(cellfun(@(x) x(trainingtrial{iRef},iB)',...
%                         spkhist(inclst{j,iclst},iRef),'UniformOutput',false));
%                     testingspk{iclst,iRef} = cell2mat(cellfun(@(x) x(testtrial{iRef},iB)',...
%                         spkhist(inclst{j,iclst},iRef),'UniformOutput',false));
%                     end
                 end
                 mdl = fitcsvm(cell2mat(trainingspk)',[ones(1,ntraining),ones(1,ntraining)*2]');
                    predicted = predict(mdl,cell2mat(testingspk)')==[ones(1,ntraining),ones(1,ntraining)*2]';
                    accuracy{iS,iclst}(iIter,iB) = nanmean(predicted);
%                 if i==1
%                     mdl = fitcsvm(cell2mat(trainingspk)',[ones(1,ntraining),ones(1,ntraining)*2]');
%                     predicted = predict(mdl,cell2mat(testingspk)')==[ones(1,ntraining),ones(1,ntraining)*2]';
%                     accuracy{iS,i}(iIter,iB) = nanmean(predicted);
%                 else
%                     mdl = fitcsvm(cell2mat(trainingspk([1:4]~=i-1,:))',[ones(1,ntraining),ones(1,ntraining)*2]');
%                     predicted = predict(mdl,cell2mat(testingspk([1:4]~=i-1,:))')==[ones(1,ntraining),ones(1,ntraining)*2]';
%                     accuracy{iS,i}(iIter,iB) = nanmean(predicted);
%                 end
            end
         end
    end    
end

ct = [0 0 0; cbrewer('qual','Dark2',3)];
time = movmean(win(1):step:win(2),5,'Endpoints','discard');
out = cellfun(@isempty,accuracy(:,1));
score = cellfun(@nanmean,accuracy(~out,:),'UniformOutput',false);
fHandle= figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 3]);
hold on;
for iClst = 1:4
    m = nanmean(cell2mat(score(:,iClst)));
    s = nanstd(cell2mat(score(:,iClst)))/sqrt(size(score,1));
   fill([time flip(time)],[m+s flip(m-s)],ct(iClst,:),'EdgeColor','none');
   plot(time,m,'Color',ct(iClst,:));
end
ylabel('Decoding accuracy');
xlabel('Time from ripple (s)');
set(gca,'Box','off','TickDir','out','FontSize',5);
alpha(0.2);
xlim([-2 2])
cd('D:\OneDrive - University of California, San Francisco\figures\allen');
print(fHandle,'-dtiff','-r600','decoding.tif');