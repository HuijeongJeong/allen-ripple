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

rippletype = {'medial';'lateral';'global'};
[spkavepreripple,delta_auc,pre_auc,post_auc,centroid_x,centroid_y] = deal(nan(nS,4,5));
spkaveconv = nan(nS,4,5,nbin);
nneuron = nan(nS,4);

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

    [spkhistz,spkpreripple] = deal(cell(1,3));    
    for iRp = 1:5
        if iRp<4
            ripple = CA1_ripple_classified.(rippletype{iRp})(:,1);
        else
            ripple = sort([CA1_ripple_classified.(rippletype{iRp-3})(:,1);...
                CA1_ripple_classified.global(:,1)]);
        end
        
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
    end
    
    spkave = cell(4,5);
    for iClst = 1:4
       spkave(iClst,:) = cellfun(@(x) nanmean(x(cidx==iClst-1,:)),spkhistz,'UniformOutput',false); 
       spkavepreripple(iS,iClst,:) = cellfun(@(x) nanmean(x(cidx==iClst-1)),spkpreripple);
       nneuron(iS,iClst) = sum(cidx==iClst-1);
    end
    if any(nneuron(iS,:)<5)
       continue; 
    end
    spkaveconv_temp = cellfun(@(x) conv2(x,fspecial('Gaussian',[1 5*resolution],resolution),'same'),spkave,'UniformOutput',false);
    pre_auc(iS,:,:) = cellfun(@(x) sum(x(time>=-1 & time<0))*binsize,spkaveconv_temp);
    post_auc(iS,:,:) = cellfun(@(x) sum(x(time>0 & time<=1))*binsize,spkaveconv_temp);
    delta_auc(iS,:,:) = abs(pre_auc(iS,:,:))-abs(post_auc(iS,:,:));
    
    [centroid_x(iS,:,:),centroid_y(iS,:,:)] = cellfun(@(x) centroid(polyshape([-1,time(time>=-1&time<=1),1],...
        [0,abs(x(time>=-1&time<=1)),0])),spkaveconv_temp);

    spkaveconv(iS,:,:,:) = permute(cell2mat(reshape(spkaveconv_temp,[4,1,5])),[1,3,2]);
end


%%
dir = 'D:\OneDrive - UCSF\figures\allen';

%%
close all
clr_l = cbrewer('qual','Pastel2',4);
clr_l = [clr_l(end,:);clr_l(1:3,:)];
clr = cbrewer('qual','Dark2',4);
clr = [clr(end,:);clr(1:3,:)];
titleList = {'dCA1';'iCA1';'global';'dCA1+global';'iCA1+global'};

data = centroid_x(:,[4,1,2,3],:);
jRp = 3;
in = sum(isnan(squeeze(data(:,:,iRp))),2)==0;
[t,p,s,m] = deal(nan(4,3));

for iRp = 1:3
    if iRp<3
        if iRp==1
        fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 5 3.5]);
        end
        axes('Position',axpt(2,1,iRp,1,axpt(10,10,2:10,1:9))); hold on;
    elseif iRp==3
        fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.9 3.5]);
        hold on;
    end
    m(:,iRp) = nanmean(squeeze(data(:,:,iRp)));
    s(:,iRp) = nanstd(squeeze(data(:,:,iRp)))/sqrt(sum(in));
    plot(1:4,squeeze(data(:,:,iRp)),'Color',[0.6 0.6 0.6]);
    errorbar(1:4,nanmean(squeeze(data(:,:,iRp))),...
        nanstd(squeeze(data(:,:,iRp)))/sqrt(sum(~isnan(data(:,1,iRp)))),'k','CapSize',3);
    plot([0,5],[0,0],'k:');
    xlim([0.5,4.5]);
    ylim([-0.3, 0.4]);
    set(gca,'XTick',1:4,'XTickLabel',{'Inh','Thal','iAct','dAct'},...
        'XTickLabelRotation',45,'YTick',-0.2:0.2:0.4,'Box','off','Tickdir','out','LineWidth',0.35);
    if rem(iRp,2)==1
        ylabel('Centroid (s)');
    else
        set(gca,'YTickLabel',[]);
    end
    if iRp==3
        x = repmat(1:4,size(data,1),1);
        y = squeeze(data(:,:,iRp));
        [r,p] = corr(x(:),y(:),'rows','complete');
%         print(fHandle,'-depsc',[dir,'\fig4\centroid_global.ai']);
    elseif iRp==2
%         print(fHandle,'-depsc',[dir,'figS1\centroid_local.ai']);
    end
end
%%
for iRp = 1:3
for iClst = 1:4
[~,p(iClst,iRp),~,stat] = ttest(squeeze(data(in,iClst,iRp)));    
t(iClst,iRp) = stat.tstat;
end
end
%%

close all
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
   [~,p,~,stat] = ttest(data(:,1),data(:,2))
   if p<0.05
       text(1.5,0.15,'*','FontSize',8);
   end
   plot([0.5 2.5],[0 0],'k:');
   ylim([-0.5 0.2]);
   xlim([0.5 2.5]);
   set(gca,'Box','off','TickDir','out','FontSize',8,'XTick',1:2,...
       'XTickLabel',{'dCA1 ripple','iCA1 ripple'},'YTick',-0.4:0.2:0.2,...
       'XTickLabelRotation',45);
   if i==2
       set(gca,'YTickLabel',[]);
   else
       ylabel({'Ripple-associated'; 'suppression'});
   end
end
print(fHandle,'-depsc',[dir,'\figS1\suppression_comparison.ai']);


%%
pre_data = squeeze(pre_auc(:,[4,1,2,3],:));
post_data = squeeze(post_auc(:,[4,1,2,3],:));
delta_data = squeeze(delta_auc(:,[4,1,2,3],:));
figure
for iRp = 1:5
   h(1) = axes('Position',axpt(5,2,iRp,1)); hold on;
   h(2) = axes('Position',axpt(5,2,iRp,2)); hold on;
   errorbar(1:4,nanmean(squeeze(pre_data(:,:,iRp))),...
       nanstd(squeeze(pre_data(:,:,iRp)))/sqrt(sum(~isnan(pre_data(:,1,iRp)))),'k','Parent',h(1));
   errorbar(1:4,nanmean(squeeze(post_data(:,:,iRp))),...
       nanstd(squeeze(post_data(:,:,iRp)))/sqrt(sum(~isnan(post_data(:,1,iRp)))),'k--','Parent',h(1));
   errorbar(1:4,nanmean(squeeze(delta_data(:,:,iRp))),...
       nanstd(squeeze(delta_data(:,:,iRp)))/sqrt(sum(~isnan(delta_data(:,1,iRp)))),'k','Parent',h(2));
   plot([0 5],[0 0],'k:','Parent',h(1));
   plot([0 5],[0 0],'k:','Parent',h(2));
   xlim([0 5]);
   ylim(h(1),[-0.4 0.4]);
   ylim(h(2),[-0.2 0.4]);
   if iRp==1
       ylabel(h(1),'AUC');
       ylabel(h(2),'\DeltaAUC(abspre-abspost)');
       legend(h(1),{'Pre';'Post'});
   end
   title(h(1),titleList{iRp});
   set(h(2),'XTick',1:4,'XTickLabel',{'Inh';'Thal';'iAct';'dAct'},'XTickLabelRotation',45);
   set(h(1),'XTick',1:4,'XTickLabel',[]);
   
end

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
