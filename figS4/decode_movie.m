rng(2)
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
meansse = nan(nS,4);

% parpool(7)
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_cellTable.mat']);
    load([sdir(sessionList(iS)),'_Events.mat'],'natural_movie');
    
    invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis & tag.celltype.rs));
%     invis = ismember(T.unit_id,tag.info.unit_id(tag.area.vis));
    sig_rf = T.p_value_rf<0.01 & T.area_rf<2500;
    inanal = invis & sig_rf;
    
    [in,idx] = ismember(T.unit_id(inanal),unitid);
    clstidx = zeros(sum(inanal),1);
    clstidx(in) = cluster_idx.vis{nclst-1}(idx(in));
    
    spktime = T.spike_time(inanal);
    spkhist = cellfun(@(y) cell2mat(cellfun(@(x) histc(y,x+[0:30])',num2cell(natural_movie.window(:,1)),...
        'UniformOutput',false)),spktime,'UniformOutput',false);
    spkhist = cellfun(@(x) x(:,1:end-1),spkhist,'UniformOutput',false);
    spkhistz = cellfun(@(x) (x-mean(x(:)))/std(x(:)),spkhist,'UniformOutpu',false);
    spkhistz = cat(3,spkhistz{:});
    
    nneurons = [sum(clstidx==0),sum(clstidx==1),sum(clstidx==2),sum(clstidx==3)];
    if min(nneurons)<5
        continue;
    end
    
    %%
    ntrial = size(natural_movie.window,1);
    sse = nan(ntrial,4);
    parfor iT = 1:ntrial
        traintrial = [1:iT-1,iT+1:ntrial];
        for iClst = 1:4
            inneuron = randsample(find(clstidx==iClst-1),min(nneurons));
            %%
            traindata = permute(squeeze(spkhistz(traintrial,:,inneuron)),[3,2,1]);
            trainidx = repmat([1:30]',1,length(traintrial));
            model = fitcecoc(traindata(:,:)',trainidx(:));
            
            %%
            testdata = squeeze(spkhistz(iT,:,inneuron));
            predictedresult = predict(model,testdata);
            %%
            sse(iT,iClst) = sqrt(sum((predictedresult-[1:30]').^2));
        end
    end
    meansse(iS,:) = mean(sse);
end

%%
close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2,2,3,4]);
plot(1:4,meansse,'Color',[0.6 0.6 0.6]);
hold on;
errorbar(1:4,nanmean(meansse),nanstd(meansse)/sqrt(sum(~isnan(meansse(:,1)))),'k','LineWidth',1);
xlim([0 5]);
ylim([20 70]);
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'YTick',20:20:60,...
    'XTick',1:4,'XTickLabel',{'Nomod','iAct','dAct','Inh'},'XTickLabelRotation',45);
ylabel('Sum of squared error');
print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\figures\2.allen\revision\FigS_visualresponse\movie_decoding_error_rs.ai');
[tbl,rm] = simple_mixed_anova(meansse); %p=0.5616



