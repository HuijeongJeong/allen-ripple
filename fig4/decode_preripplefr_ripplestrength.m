clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

basewin = [-3 -2];
win = [-1 0];
rippletype = {'medial';'lateral';'global'};

winhist = [-3 3];
binsize = 0.01;
resolution = 5;

[preripplespk_strength,preripplespk_length] = deal(nan(nS,3,5,3));
[r_strength,p_strength,r_length,p_length] = deal(nan(nS,3,5));
[averipple_strength,averipple_length] = deal(nan(nS,3,3));
[corr_strength,corr_length] = deal(nan(nS,3,10,3));
avespkhistz = nan(nS,3,5,3,diff(winhist)/binsize);

nIter = 100;
nCell = 100;
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_ripple_modulated.mat'],'fr_ripple');
    load([sdir(sessionList(iS)),'_cellTable.mat'],'T');
    load([sdir(sessionList(iS)),'_ripples.mat'],'CA1_ripple_classified');
    
    targets = tag.info.unit_id((tag.area.vis&tag.celltype.rs) | tag.area.thalamus);
    in = ismember(T.unit_id,targets);
    
    unitid = T.unit_id(in);
    spiketime = T.spike_time(in);
    
    cidx = zeros(sum(in),1);
    in = ismember(unitid,tag.info.unit_id(tag.area.thalamus));
    cidx(in) = 4; % thalamus
    [in,idx] = ismember(unitid,unit_id.vis);
    cidx(in) = cluster_idx.vis{2}(idx(in));
    
    clstnum = cellfun(@(x) sum(cidx==x),num2cell(0:4));
    if min(clstnum)<5
        continue;
    end
    
    for iRp = 1:3
        preripplespk = cellfun(@(y) cellfun(@(x) sum(y>=x+win(1) & y<=x+win(2)),...
            num2cell(CA1_ripple_classified.(rippletype{iRp})(:,1)))',...
            spiketime,'UniformOutput',false);
        preripplespkz = cell2mat(cellfun(@(x) zscore(x,[],2),preripplespk,'UniformOutput',false));
%         baseripplespk = cellfun(@(y) cellfun(@(x) sum(y>=x+basewin(1) & y<=x+basewin(2)),...
%             num2cell(CA1_ripple_classified.(rippletype{iRp})(:,1)))',...
%             spiketime,'UniformOutput',false);
%         preripplespkz = cell2mat(cellfun(@(x,y) (x-mean(y))/std(y),...
%             preripplespk,baseripplespk,'UniformOutput',false));
%         preripplespkz(isinf(preripplespkz)) = nan;
        
        if iRp==3
            ripplestrength = max(CA1_ripple_classified.global_strength,[],2);
        else
            ripplestrength = CA1_ripple_classified.(rippletype{iRp})(:,4);
        end
        
        nRipple = length(ripplestrength);
        predictedripple = nan(nIter,6);
        actualripple = nan(nIter,1);
        betacells = nan(nIter,5);
        for iIter = 1:nIter
            iIter
            testingripple = randsample(nRipple,1);
            trainingripples = [1:testingripple-1,testingripple+1:nRipple];
            actualripple(iIter) = ripplestrength(testingripple);
            for iClst = 1:5
                trainingcells = find(cidx==iClst-1);
                %                 trainingcells = randsample(find(cidx==iClst-1),min(clstnum));
                %                 mdl = fitrsvm(preripplespkz(trainingcells,trainingripples)',...
                %                     ripplestrength(trainingripples));
                net = fitnet(10,'trainlm');
                net.divideParam.trainRatio = 75/100;
                net.divideParam.valRatio = 15/100;
                net.divideParam.testRatio = 15/100;
                [net,tr] = train(net,preripplespkz(trainingcells,:),ripplestrength');
                y = net(preripplespkz(trainingcells,:));
                rss = sum((y(tr.testInd)-ripplestrength(tr.testInd)').^2);
                tss = sum((ripplestrength(tr.testInd)-mean(ripplestrength(tr.testInd))).^2);
                rsquare(iIter,iClst) = 1-rss/tss;
                
                mdl = fitglm(preripplespkz(trainingcells,trainingripples)',...
                    ripplestrength(trainingripples));
                predictedripple(iIter,iClst) = predict(mdl,preripplespkz(trainingcells,testingripple)');
            end
            
            trainingcells = 1:length(unitid);
%             mdl = fitrsvm(preripplespkz(trainingcells,trainingripples)',...
%                 ripplestrength(trainingripples));
%             mdl = fitglm(preripplespkz(trainingcells,trainingripples)',...
%                 ripplestrength(trainingripples));
            net = fitnet(10,'trainlm');
            net.divideParam.trainRatio = 75/100;
            net.divideParam.valRatio = 15/100;
            net.divideParam.testRatio = 15/100;
            [net,tr] = train(net,preripplespkz(trainingcells,:),ripplestrength');
            y = net(preripplespkz(trainingcells,:));
            rss = sum((y(tr.testInd)-ripplestrength(tr.testInd)').^2);
            tss = sum((ripplestrength(tr.testInd)-mean(ripplestrength(tr.testInd))).^2);
            rsquare(iIter,6) = 1-rss/tss;
%             predictedripple(iIter,end) = predict(mdl,preripplespkz(trainingcells,testingripple)');
%             for iClst = 1:5
%                inclst = cidx(trainingcells)==iClst-1;
%                betacells(iIter,iClst) = mean(mdl.Beta(inclst));
%             end
        end
    end
end

%% [0,1],[0,2],[0,3],[0,4],[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]
%%
dir = 'D:\OneDrive - University of California, San Francisco\figures\allen';

for iRp = 1:3
    if iRp<3
        if iRp==1
            fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5 3.5]);
        end
        subplot(1,2,iRp)
    elseif iRp==3
        fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3.75 3.5]);
    end
    hold on;
    plot(1:4,squeeze(r_strength(:,iRp,[5,2,3,4])),'color',[0.6 0.6 0.6],'LineWidth',0.35);
    errorbar(1:4,nanmean(squeeze(r_strength(:,iRp,[5,2,3,4]))),...
        nanstd(squeeze(r_strength(:,iRp,[5,2,3,4])))/sqrt(sum(~isnan(r_strength(:,1,1)))),...
        'k','CapSize',3);
    [h,pttest(iRp,:)] = ttest(squeeze(r_strength(:,iRp,[5,2,3,4])));
    for iclst = 1:4
       if h(iclst)==1
           text(iclst,0.5,'*','FontSize',8);
       end
    end
    plot([0 6],[0 0],'k:');
    ylim([-0.6 0.6])
    xlim([0.5 4.5])
    set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',...
        1:4,'XTickLabel',{'Thal','iAct','dAct','Inh'},'XTickLabelRotation',45,'YTick',-0.6:0.3:0.6);
    if rem(iRp,2)==1
        ylabel({'Correlation'; '(pre-ripple firing'; 'vs. ripple strength)'});
    else
        set(gca,'YTickLabel',[]);
        print(fHandle,'-depsc',[dir,'\corr_prespk_ripplestrength_localizedripples.ai']);
    end
end
print(fHandle,'-depsc',[dir,'\fig4\corr_prespk_ripplestrength_globalripples.ai']);

%%
jClst = [5,2,4];
time = winhist(1)+binsize/2:binsize:winhist(2)-binsize/2;
ct = [cbrewer('qual','Set2',4);cbrewer('qual','Dark2',4)];
for iRp = 1:3
    if iRp<3
        if iRp==1
            fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5 3.5]);
        end
    elseif iRp==3
        fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3.75 3.5]);
    end
    for iClst = 1:3
        if iRp<3
        axes('Position',axpt(3,2,iClst,iRp)); hold on;
        else
            axes('Position',axpt(3,1,iClst,1)); hold on;
        end
        for i = 1:3
        data = squeeze(avespkhistz(:,iRp,jClst(iClst),i,:));
        plot(time,nanmean(data));   
        end
    end
end

%%
close all

jClst = [5,2,3,4];
ct = [cbrewer('qual','Set2',4);cbrewer('qual','Dark2',4)];

for iRp = 1:3
    if iRp<3
        if iRp==1
            fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 7 8]);
        end
    elseif iRp==3
        print(fHandle,'-depsc',[dir,'\ave_prespk_ripplestrength_localizedripples.ai']);
        fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 7 3.5]);
    end
    for iClst = 1:4
        if iRp<3
            axes('Position',axpt(4,2,iClst,iRp,axpt(10,10,3:10,1:8))); hold on;
        else
            axes('Position',axpt(4,1,iClst,1,axpt(10,10,3:10,1:8))); hold on;
        end
        data = squeeze(preripplespk_strength(:,iRp,jClst(iClst),:));
        plot(1:3,data,'Color',ct(jClst(iClst)-1,:),'LineWidth',0.35);
        errorbar(1:3,nanmean(data),nanstd(data)/sqrt(sum(~isnan(data(:,1)))),...
            'Color',ct(jClst(iClst)-1+4,:),'LineWidth',0.75,'CapSize',3);
        plot([0.5 3.5],[0 0],'k:','LineWidth',0.35);
        xlim([0.5 3.5]);
        ylim([-1 1]);
        set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
            'XTick',1:3,'YTick',-1:0.5:1,'XTickLabel',{'Bottom 25%';'Middle';'Top 25%'},...
            'XTickLabelRotation',45);
        if iClst==1
           ylabel({'Norm. firing rate'; '(z-score)'});
           if iRp~=1
           xlabel('Ripple power');
           end
        else
           set(gca,'YTickLabel',[]); 
        end
        if iRp==1
           set(gca,'XTickLabel',[]); 
        end
    end
end
print(fHandle,'-depsc',[dir,'\fig4\ave_prespk_ripplestrength_globalripples.ai']);
% print(fHandle,'-depsc',[dir,'\fig4\corr_prespk_ripplestrength_globalripples.ai']);

