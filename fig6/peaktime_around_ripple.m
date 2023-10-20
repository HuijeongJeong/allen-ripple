clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

winbase = [-3, -2];
wintest = [-1, 0];
win = [-3 3];
binsize = 0.01;
time = win(1)+binsize/2:binsize:win(2)-binsize/2;
resolution = 5;
nbin = diff(win)/binsize;

data = cell(nS,3);
data_g = cell(nS,2);

rippletype = {'medial';'lateral';'global'};
searchpeak = [0,1,1,0];
[maxlatency,maxvalue,m50latency,spkavepreripple] = deal(nan(nS,4,3));
spkaveconv = nan(nS,4,3,nbin);
nneuron = nan(nS,4);
nssave = nan(nS,2,3,nbin);

for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_ripple_modulated.mat'],'fr_ripple');
    load([sdir(sessionList(iS)),'_cellTable.mat'],'T');
    load([sdir(sessionList(iS)),'_ripples.mat'],'CA1_ripple_classified','filtered_lfp');
    
    targets = tag.info.unit_id((ismember(tag.info.unit_id,unit_id.vis)&tag.celltype.rs) |...
        tag.area.thalamus);
    in = ismember(T.unit_id,targets);
    
    unitid = T.unit_id(in);
    spiketime = T.spike_time(in);
    
    cidx = zeros(sum(in),1);
    [in,idx] = ismember(unitid,unit_id.vis);
    cidx(in) = cluster_idx.vis{2}(idx(in));

    [nss,nsstime] = deal(cell(2,1));
    for iProbe = 1:2
        if iscell(filtered_lfp.time{iProbe})
            nss_tmp = cell(length(filtered_lfp.time{iProbe}),1);
            for iw = 1:length(filtered_lfp.time{iProbe})
                nss_tmp{iw} =  mean(NormalizedSquaredSignal_HJ([filtered_lfp.time{iProbe}{iw}',...
                    filtered_lfp.lfp{iProbe}{iw}]),2);
            end
            nss{iProbe} = cell2mat(nss_tmp);
            nsstime{iProbe} = cell2mat(filtered_lfp.time{iProbe}');
        else
            nss{iProbe} = mean(NormalizedSquaredSignal_HJ([filtered_lfp.time{iProbe}',...
                filtered_lfp.lfp{iProbe}]),2);
            nsstime{iProbe} = filtered_lfp.time{iProbe};
        end
    end
        
    [spkhistz,spkpreripple] = deal(cell(1,3));    
    rippletotal = sort([CA1_ripple_classified.medial(:,1);CA1_ripple_classified.lateral(:,1);...
        CA1_ripple_classified.global(:,1)]);
    for iRp = 1:3
        
%         outripple = cellfun(@(x) x-rippletotal(find(x-rippletotal>0,1,'last')),...
%             num2cell(CA1_ripple_classified.(rippletype{iRp})(:,1)),'UniformOutput',false);
%         outripple(cellfun(@isempty,outripple)) = {nan};
%         outripple = cell2mat(outripple)<1;
%         ripple = CA1_ripple_classified.(rippletype{iRp})(~outripple,1);
        ripple = CA1_ripple_classified.(rippletype{iRp})(:,1);
        
        spkbaseline = cellfun(@(y) cellfun(@(x) sum(y>=x+winbase(1) & y<=x+winbase(2)),...
            num2cell(ripple))',spiketime,'UniformOutput',false);
        spktest = cellfun(@(y) cellfun(@(x) sum(y>=x+wintest(1) & y<=x+wintest(2)),...
            num2cell(ripple))',spiketime,'UniformOutput',false);

        
        spkhist_temp = cellfun(@(y) nanmean(cell2mat(cellfun(@(x) histc(y,x+[win(1):binsize:win(2)])'*(1/binsize),...
            num2cell(ripple),'UniformOutput',false))),spiketime,'UniformOutput',false);
        spkhistz{iRp} = cell2mat(cellfun(@(x,y) (x(1:end-1)-mean(y))/std(y),spkhist_temp,spkbaseline,'UniformOutput',false));
        spkhistz{iRp}(isinf(spkhistz{iRp})) = nan; 

        spkpreripple{iRp} = cellfun(@(x,y) (mean(x/diff(wintest))-mean(y/diff(winbase)))/...
            std(y/diff(winbase)),spktest,spkbaseline);
        spkpreripple{iRp}(isinf(spkpreripple{iRp})) = nan;
        
        for iProbe = 1:2
            intime = cellfun(@(x) nsstime{iProbe}>=x+win(1)-0.1 &...
                nsstime{iProbe}<=x+win(2)+0.1,...
                num2cell(ripple),'UniformOutput',false);
            nssave(iS,iProbe,iRp,:) = mean(cell2mat(cellfun(@(x,y) interp1(nsstime{iProbe}(y),...
                nss{iProbe}(y),x+[win(1)+binsize/2:binsize:win(2)-binsize/2]),...
                num2cell(ripple),intime,'UniformOutput',false)));
        end
    end
    
    spkave = cell(4,3);
    for iClst = 1:4
       spkave(iClst,:) = cellfun(@(x) nanmean(x(cidx==iClst-1,:)),spkhistz,'UniformOutput',false); 
       spkavepreripple(iS,iClst,:) = cellfun(@(x) nanmean(x(cidx==iClst-1)),spkpreripple);
       nneuron(iS,iClst) = sum(cidx==iClst-1);
    end
    if any(nneuron(iS,:)<5)
       continue; 
    end
    spkaveconv_temp = cellfun(@(x) conv2(x,fspecial('Gaussian',[1 5*resolution],resolution),'same'),spkave,'UniformOutput',false);
    
    latencyidx = nan(4,3);
    [maxvalue(iS,searchpeak==1,:),latencyidx(searchpeak==1,:)] = cellfun(@(x) max(x),spkaveconv_temp(searchpeak==1,:));
    [maxvalue(iS,searchpeak==0,:),latencyidx(searchpeak==0,:)] = cellfun(@(x) min(x),spkaveconv_temp(searchpeak==0,:));
    maxlatency(iS,:,:) = time(latencyidx);
    [~,latencyidx(searchpeak==1,:)] = cellfun(@(x,y) find(x>y/2,1,'first'),spkaveconv_temp(searchpeak==1,:),num2cell(squeeze(maxvalue(iS,searchpeak==1,:))));
    [~,latencyidx(searchpeak==0,:)] = cellfun(@(x,y) find(x<y/2,1,'first'),spkaveconv_temp(searchpeak==0,:),num2cell(squeeze(maxvalue(iS,searchpeak==0,:))));
    m50latency(iS,:,:) = time(latencyidx);    
    
    spkaveconv(iS,:,:,:) = permute(cell2mat(reshape(spkaveconv_temp,[4,1,3])),[1,3,2]);
end

%%
close all
dir = 'D:\OneDrive - University of California, San Francisco\figures\allen';
auc = sum(spkaveconv(:,:,:,time>=-0.5&time<0.5),4)*binsize;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 3.5]);
for i = 1:2
   axes('Position',axpt(2,1,i,1,axpt(10,10,2:9,1:9),[0.05 0.05]))
   hold on;
   if i==1
       data = squeeze(auc(:,1,1:2));
       title('Thalamus');
   else
       data = squeeze(auc(:,4,1:2));
       title('Inh');
   end
   plot(1:2,data,'Color',[0.6 0.6 0.6]);
   errorbar(1:2,nanmean(data),nanstd(data)/sqrt(sum(~isnan(data(:,1)))),'k');
   [~,p] = ttest(data(:,1),data(:,2));
   if p<0.05
       text(1.5,0.15,'*','FontSize',8);
   end
   plot([0.5 2.5],[0 0],'k:');
   ylim([-0.5 0.2]);
   xlim([0.5 2.5]);
   set(gca,'Box','off','TickDir','out','FontSize',8,'XTick',1:2,...
       'XTickLabel',{'dCA1 ripple','iCA1 ripple'},'YTick',-0.4:0.2:0.2,...
       'XTickLabelRotation',45,'LineWidth',0.35);
   if i==2
       set(gca,'YTickLabel',[]);
   else
       ylabel({'Ripple-associated'; 'suppression'});
   end
end
print(fHandle,'-depsc',[dir,'\figS1\suppression_comparison.ai']);

%%
close all; figure;
ct = [cbrewer('qual','Set2',4);cbrewer('qual','Dark2',4)];
x = [4,1,2,3];
rippletype = {'dCA1','iCA1','Global'};
dir = 'D:\OneDrive - University of California, San Francisco\figures\allen';
for iRp = 1:3
    if iRp<3
        if iRp==1
        fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 8 3.5]);
        end
        axes('Position',axpt(2,1,iRp,1,axpt(10,10,2:9,1:9),[0.1 0.05]))
    elseif iRp==3
        fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5 3.5]);
    end
    hold on;
    
    colororder({'k','c'})
    yyaxis right
    if iRp==3
        ylabel({'Norm. ripple';'band power'})
        m = nanmean(squeeze(mean(nssave(:,:,iRp,:),2)),1);
        plot(time,m,'LineWidth',1);
        set(gca,'YTick',0:5:10,'YLim',[-2 10]);
    else
        yyaxis right
        plot(win(1)+binsize/2:binsize:win(2)-binsize/2,nanmean(squeeze(nssave(:,iRp,iRp,:)),1));
        if iRp==2
            set(gca,'YTick',0:5:10,'YLim',[-2 12]);
            ylabel({'Norm. ripple';'band power'})
        else
            set(gca,'YTick',0:3:6,'YLim',[-1 7]);
        end
    end
    
    yyaxis left
    for iClst = 1:4
        m = nanmean(squeeze(spkaveconv(:,iClst,iRp,:)),1);
        s = nanstd(squeeze(spkaveconv(:,iClst,iRp,:)),1)/sqrt(sum(~isnan(spkaveconv(:,1,1,1))));
        %         plot(time,squeeze(spkaveconv(:,iClst,iRp,:)),'Color',ct(x(iClst),:),'LineWidth',0.35);
        fill([time,flip(time)],[m+s flip(m-s)],ct(x(iClst),:),'EdgeColor','none');
        plot(time,m,'Color',ct(x(iClst)+4,:),'LineWidth',1,'LineStyle','-','Marker','none');
    end
    plot([0 0],[-0.7 1.2],'k:');
    ylim([-0.5 0.8])
    xlim([-2.5 2.5])   
    set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',...
        -2:2:2,'YTick',-0.4:0.4:0.8);
    xlabel('Time (s)');
    if rem(iRp,2)==1
        ylabel({'Norm. firing rate';'(z-score)'});
    else
            set(gca,'YTickLabel',[]);
%             print(fHandle,'-depsc',[dir,'\avepsth_localizedripples.ai']);
    end
end
% print(fHandle,'-depsc',[dir,'\fig4\avepsth_globalripples.ai']);

%%
close all; 
for iRp = 1:3
    if iRp<3
        if iRp==1
        fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5 3.5]);
        end
        subplot(1,2,iRp)
    elseif iRp==3
        fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 3 3.5]);
    end
    hold on;
    data = squeeze(m50latency(:,:,iRp));
    plot(1:4,data,'Color',[0.6 0.6 0.6],'LineWidth',0.35);
    errorbar(1:4,nanmean(data),nanstd(data)/sqrt(sum(~isnan(data(:,1)))),...
        'k','LineWidth',0.5,'CapSize',3);
    plot([0 5],[0 0],'k:')
    ylim([-3 0.5]);
    xlim([0.5 4.5])
    set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,'XTick',...
        1:4,'YTick',-3:1,'XTickLabel',{'Thal';'iAct';'dAct';'Inh'},'XTickLabelRotation',45);
    if rem(iRp,2)==1
       ylabel({'Latency from'; 'ripple onset (s)'});
    else
        set(gca,'YTickLabel',[]);
%         print(fHandle,'-depsc',[dir,'\latency_50_localizedripples.ai']);
    end
end
% % print(fHandle,'-depsc',[dir,'\fig4\latency_50_globalripples.ai']);
