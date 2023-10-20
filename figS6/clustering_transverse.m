rng(2)

clearvars; clc; close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric','units');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.global_ripple_number>0 & session_metric.distance_bw_ml_probes>2000);
nS = length(sessionList);

overlap_fraction = nan(nS,5,2);
data = cell(nS,5,2);
unitid = cell(nS,1);
distance2medial = cell(nS,5);

win = [-5,5];
analwin = [-1.5,1.5];
bin = 0.01;
resolution = 10;
time = win(1)+bin/2:bin:win(2)-bin/2;
intime = time>=analwin(1) & time<=analwin(2);

%%
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_win',...
        'spontaneous_anal_win','spontaneous_CA1_ripple','CA1_ripple_classified');
    load([sdir(sessionList(iS)),'_cellTable.mat'],'T');
    load([sdir(sessionList(iS)),'_ripple_modulated.mat'],'fr_ripple');
    
 
    n = size(spontaneous_CA1_ripple,1); 
    medialidx = find(strcmp(spontaneous_CA1_ripple.relative_location,'medial'));
    lateralidx = find(strcmp(spontaneous_CA1_ripple.relative_location,'lateral'));
    ccf = cell2mat(spontaneous_CA1_ripple.ccf_coordinate);
    d2m = abs(ccf-ccf(medialidx,:));
    [~,sortidx] = sort(sum(d2m,2));    
    distance2medial(iS,1:n-1) = mat2cell(d2m(sortidx(2:end),:),ones(n-1,1),3);

    ripples = spontaneous_CA1_ripple.ripple(sortidx);
    ripple_ref = ripples{1};

%%
    inunit = ismember(T.unit_id,unit_id.vis);
    unitid{iS} = T.unit_id(inunit);
    
    %%
    for irp = 1:n-1
        ripple = cell(1,2);
        ripple_compare = ripples{irp+1};
        overlapidx_m2d = cellfun(@(x,y) sum((x>=ripple_compare(:,1) & x<=ripple_compare(:,3)) |...
            (y>=ripple_compare(:,1) & y<=ripple_compare(:,3)) |...
            (x<=ripple_compare(:,1) & y>=ripple_compare(:,3)))>0,...
            num2cell(ripple_ref(:,1)),num2cell(ripple_ref(:,3)));
        overlapidx_d2m = cellfun(@(x,y) sum((x>=ripple_ref(:,1) & x<=ripple_ref(:,3)) |...
            (y>=ripple_ref(:,1) & y<=ripple_ref(:,3)) |...
            (x<=ripple_ref(:,1) & y>=ripple_ref(:,3)))>0,...
            num2cell(ripple_compare(:,1)),num2cell(ripple_compare(:,3)));
        ripple{1} = ripple_ref(~overlapidx_m2d,:);
        ripple{2} = ripple_compare(~overlapidx_d2m,:);
        
        overlapping_ripples_m = ripple_ref(overlapidx_m2d,:);
        overlapping_ripples_l = ripple_compare(overlapidx_d2m,:);
        for i = 1:length(overlapping_ripples_m)
            x = overlapping_ripples_l(:,1);
            y = overlapping_ripples_l(:,3);
            idx = find((x>=overlapping_ripples_m(i,1) & x<=overlapping_ripples_m(i,3))|...
                (y>=overlapping_ripples_m(i,1) & y<=overlapping_ripples_m(i,3)) |...
                (x<=overlapping_ripples_m(i,1) & y>=overlapping_ripples_m(i,3))>0);
            if overlapping_ripples_m(i,4)>max(overlapping_ripples_l(idx,4))
                overlapping_ripples_m(i,5) = 1;
            else
                overlapping_ripples_l(idx,5) = 1;
            end
        end
%         ripple{3} = [ripple{1};overlapping_ripples_m(overlapping_ripples_m(:,5)==1,1:4)];
%         ripple{4} = [ripple{2};overlapping_ripples_l(overlapping_ripples_l(:,5)==1,1:4)];
              
        spktime = cellfun(@(y) cellfun(@(x) spikeWin(x,y(:,1),win),...
            T.spike_time(inunit),'UniformOutput',false),ripple,'UniformOutpu',false);
        spkhist = cellfun(@(z) cellfun(@(y) cell2mat(cellfun(@(x) histcounts(x,win(1):bin:win(2))*(1/bin),...
            y,'UniformOutput',false)),z,'UniformOutput',false),spktime,'UniformOutput',false);
        spkconv = cellfun(@(z) cell2mat(cellfun(@(x) conv(mean(x,1),...
            fspecial('Gaussian',[1, 5*resolution],resolution),'same'),z,...
            'UniformOutput',false)),spkhist,'UniformOutput',false);
        
        overlap_fraction(iS,irp,:) = [mean(overlapidx_m2d),mean(overlapidx_d2m)];
        data(iS,irp,:) = cellfun(@(x) x(:,intime),spkconv,'UniformOutput',false);
     end
end

%%

gap = nan(5,5);
avesil = nan(5,5);
f = nan(5,6);
probMat = cell(5,6);
expvar = nan(5,6);
ncell = nan(5,1);

%%
nIter = 100;
nprobes = sum(cellfun(@(x) ~isempty(x),squeeze(data(:,:,1))),2);
in = nprobes>=4;


%%
close all
for iref = 1:2
    data_sub = [];
for irp = 1:4
    if irp<3
        data_sub = [data_sub,overlap_fraction(in,irp,iref)]; 
    else
        temp = cellfun(@(x,y) squeeze(x(:,y-(4-irp),:))',mat2cell(overlap_fraction(in,:,iref),...
            ones(sum(in),1),5),num2cell(nprobes(in)),'UniformOutput',false);
        data_sub = [data_sub,cat(1,temp{:})]; 
    end
end
subplot(1,2,iref);
if iref==1
    x = 1:4;
else
    x = 4:-1:1;
end    
plot(x,data_sub,'Color',[0.6 0.6 0.6]);
hold on;
errorbar(x,mean(data_sub),std(data_sub)/sqrt(sum(in)),'k','CapSize',3);
ylim([0,1])
xlim([0,5]);
set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
    'XTick',1:4,'YTick',0:0.2:1,'XTickLabel',{'dCA1','','','iCA1'},...
    'XTickLabelRotation',45);
xlabel('Location of other probe');
if iref==1
    ylabel('P(global|ripples at dCA1)');
else
    ylabel('P(global|ripples at iCA1)');
end    
end



%%
for irp = 1:5
    if irp<3
        data_sub = data(in,irp,1:2);
    elseif irp==5
        data_sub = cellfun(@(x,y) squeeze(x(:,y,:))',mat2cell(data(:,:,1:2),ones(nS,1),5,2),...
            num2cell(nprobes),'UniformOutput',false);
        data_sub = cat(1,data_sub{:});
    else
        data_sub = cellfun(@(x,y) squeeze(x(:,y-(4-irp),:))',mat2cell(data(in,:,1:2),ones(sum(in),1),5,2),...
            num2cell(nprobes(in)),'UniformOutput',false);
        data_sub = cat(1,data_sub{:});
    end
    
    datatotal = cellfun(@(x,y) [x,y],data_sub(:,1),data_sub(:,2),'UniformOutput',false);
    datatotal = cat(1,datatotal{:});
    dataz = zscore(datatotal,[],2);
    [~,score,~,~,explained] = pca(dataz);
    nPC = find(cumsum(explained)>80,1,'first');
    ncell(irp) = size(dataz,1);
    %%
    clusterIdx = nan(ncell(irp),6);
    [w,wnull,s] = deal(nan(6,1));
    tot_withinss = nan(nIter,6);
    
    kmtmp = [];
    for k = 2:6
        for iIter = 1:nIter
            [kmtmp(k-1,:,iIter),~,~,D] = kmeans(score(:,1:nPC),k);
            dist = min(D,[],2);
            tot_withinss(iIter,k) = sum(dist);
        end
        
        probmat = NaN(ncell(irp),ncell(irp));
        for iC = 1:ncell(irp)
            for jC = iC:ncell(irp)
                probmat(iC,jC) = sum(kmtmp(k-1,iC,:)==kmtmp(k-1,jC,:));
                probmat(jC,iC) = sum(kmtmp(k-1,iC,:)==kmtmp(k-1,jC,:));
            end
        end
        
        Z = linkage(probmat,'average','chebychev');
        km = cluster(Z,'MaxClust',k);
        [~,~,outperm] = dendrogram(Z,0);
        %%
        probMat{irp,k} = probmat(outperm,outperm);
        probIdx = km(outperm);
        clusterIdx(:,k) = km;
        
        avesil(irp,k) = mean(silhouette(score(:,1:nPC),clusterIdx(:,k)));
        [w(k),wnull(k),s(k)] = gapstat(score(:,1:nPC),clusterIdx(:,k),k,10);
    end
    if irp==1
        clusterIdx_1 = clusterIdx(:,2);
    end
    
    totss = sum(pdist2(score(:,1:nPC),nanmean(score(:,1:nPC))).^2);
    expvar(irp,:) = ((totss-mean(tot_withinss,1))/totss)*100;
    %%
    gap(irp,:) = wnull(2:6)-w(2:6);
    clusterIdx(:,1) = 1;
%     f(irp,:) = fk(score(:,1:nPC),clusterIdx,6);
end
%%
coclusteringIdx = cellfun(@(x,y) sum(abs((x(:)/100)-0.5))/((y^2)/2),...
    probMat,num2cell(repmat(ncell,1,6)));
coclusteringIdx(:,1) = [];
f(:,1) = [];
ct = [cbrewer('div','RdBu',9);0,0,0];
ct(3:7,:) = [];
%%
close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 9 3]);
axes('Position',axpt(4,10,3,1:9,[],[0.1 0.05]));
hold on;
for irp = 1:4
    plot(2:6,coclusteringIdx(irp,:),'Color',ct(irp,:));
    [~,maxidx] = max(coclusteringIdx(irp,:));
    scatter(maxidx+1,coclusteringIdx(irp,maxidx),4,ct(irp,:),'filled');
end
ylabel('Coclustering index');
xlim([1 7])
ylim([0.78 1.02])
xlabel('k');
set(gca,'Box','off','TickDir','out','FontSize',8,'XTick',2:6,'YTick',0.8:0.1:1);


axes('Position',axpt(4,10,4,1:9,[],[0.1 0.05]));
hold on;
for irp = 1:4
    if irp<3
        dis_sub = cell2mat(distance2medial(in,irp));
    else
        dis_sub = cell2mat(cellfun(@(x,y) x{y-(4-irp)},...
            mat2cell(distance2medial(in,:),ones(sum(in),1),5),...
            num2cell(nprobes(in)),'UniformOutput',false));
    end
    dis_sub = mean(dis_sub);
    dis_sub = sqrt(sum(dis_sub.^2));
    
    [~,maxidx] = max(coclusteringIdx(irp,:));
    scatter(dis_sub/1000,avesil(irp,maxidx+1),4,ct(irp,:));
end
ylabel('Silhouette value at best k');
% xlim([0.5 4.5])
set(gca,'Box','off','TickDir','out','FontSize',8,'YTick',0.2:0.1:0.4,...
    'YLim',[0.2 0.4],'XLim',[0 3],'XTick',0:3);
xlabel('Distance (mm)');

axes('Position',axpt(4,10,1,1:9,[],[0.1 0.05]));
hold on;
for irp = 1:4
    if irp<3
        dis_sub = cell2mat(distance2medial(in,irp));
    else
        dis_sub = cell2mat(cellfun(@(x,y) x{y-(4-irp)},...
            mat2cell(distance2medial(in,:),ones(sum(in),1),5),...
            num2cell(nprobes(in)),'UniformOutput',false));
    end
    errorbar(1:3,mean(dis_sub(:,[1,3,2])/1000),std(dis_sub(:,[1,3,2])/1000)/sqrt(sum(in)),...
        'Color',ct(irp,:),'CapSize',3);
end
xlim([0 4])
ylabel('Distance (mm)');
set(gca,'Box','off','TickDir','out','FontSize',8,'XTick',1:3,...
    'XTickLabel',{'AP','ML','DV'},'XTickLabelRotation',45,'YTick',0:1:2);

axes('Position',axpt(4,10,2,1:9,[],[0.1 0.05]));
hold on;
data_sub = [];
dis_sub = [];
for irp = 1:4
    if irp<3
        dis_sub = [dis_sub,sqrt(sum(cell2mat(distance2medial(in,irp)).^2,2))];
        data_sub = [data_sub,overlap_fraction(in,irp,1)];
        
    else
        dis_sub = [dis_sub,sqrt(sum(cell2mat(cellfun(@(x,y) x{y-(4-irp)},...
            mat2cell(distance2medial(in,:),ones(sum(in),1),5),...
            num2cell(nprobes(in)),'UniformOutput',false)).^2,2))];
        temp = cellfun(@(x,y) squeeze(x(:,y-(4-irp),:))',mat2cell(overlap_fraction(in,:,1),...
            ones(sum(in),1),5),num2cell(nprobes(in)),'UniformOutput',false);
        data_sub = [data_sub,cat(1,temp{:})];
    end
end
%%
for is = 1:sum(in)
plot(dis_sub(is,:)/1000,data_sub(is,:),'Color',[0.6 0.6 0.6]);
end
plot(mean(dis_sub)/1000,mean(data_sub),'k');
%%
for irp = 1:4
   errorbar(mean(dis_sub(:,irp))/1000,mean(data_sub(:,irp)),std(data_sub(:,irp))/sqrt(sum(in)),...
       'Color',ct(irp,:),'CapSize',3);
   errorbar(mean(dis_sub(:,irp))/1000,mean(data_sub(:,irp)),std(dis_sub(:,irp)/1000)/sqrt(sum(in)),...
       'horizontal','Color',ct(irp,:),'CapSize',3);
end
xlim([-0.2 3.3])
ylim([0 1])
ylabel('P(global|ripples at dCA1)');
xlabel('Distance (mm)');
set(gca,'Box','off','TickDir','out','FontSize',8,'XTick',0:3,'YTick',0:0.5:1);
% print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\figures\2.allen\revision\FigS_transverse\clustering_idx_distance_2.ai')
%%
close all
irp = 1;
data_sub = data(in,irp,:);
datatotal = cellfun(@(x,y) [x,y],data_sub(:,1),data_sub(:,2),'UniformOutput',false);
datatotal = cat(1,datatotal{:});
dataz = zscore(datatotal,[],2);

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 3]);
for i = 1:2
   axes('Position',axpt(2,1,i,1,axpt(10,10,2:10,1:9)));
   hold on;
   m = mean(dataz(clusterIdx_1==i,:));
   s = std(dataz(clusterIdx_1==i,:))/sqrt(sum(clusterIdx_1==i));
   fill([time(intime), flip(time(intime))],[m(1:300)+s(1:300), flip(m(1:300)-s(1:300))],...
       [1 0.6 0.6],'EdgeColor','none');
   fill([time(intime), flip(time(intime))],[m(301:600)+s(301:600), flip(m(301:600)-s(301:600))],...
       [0.6 0.6 0.6],'EdgeColor','none');
   plot(time(intime),m(1:300),'r');
   plot(time(intime),m(301:600),'k');
   plot([0 0],[-0.75 1.5],'k:');
   title(['Cluster ',num2str(i)]);
   set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
       'YTick',-0.5:0.5:1.5,'YLim',[-0.75 1.5],'XLim',[-1.5 1.5],'XTick',-1:1);
   if i==1
       ylabel({'Norm. firing rate';'(z-score)'});
   else
       set(gca,'YTickLabel',[]);
   end
   xlabel('Time from ripple onset (s)');
end
print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\figures\2.allen\revision\FigS10\clustering_dorsal_probes.ai')

%%
figure
dataz_ave = mean(reshape(dataz,[size(dataz,1),300,2]),3);
[~,maxidx] = max(dataz_ave(:,100:200),[],2);
[~,minidx] = min(dataz_ave(:,100:200),[],2);
peakidx = nan(size(clusterIdx_1));
peakidx(clusterIdx_1==1) = maxidx(clusterIdx_1==1);
peakidx(clusterIdx_1==2) = minidx(clusterIdx_1==2);
[~,sortidx] = sortrows([clusterIdx_1,peakidx]);
for i = 1:2
axes('Position',axpt(2,1,i,1));
imagesc(dataz(sortidx,[1:300]+300*(i-1)));
end
%%

function [w,wnull,s] = gapstat(data,cluster,k,niter)
n = size(data,1);
minvals = min(data);
maxvals = max(data);

withinnss_null = nan(niter,1);
for iiter = 1:niter
    randdata = [];
    for id = 1:size(data,2)
        randdata = [randdata, minvals(id)+(maxvals(id)-minvals(id))*rand(n,1)];
    end
    [~,~,~,D] = kmeans(randdata,k);
    dist = min(D,[],2);
    withinnss_null(iiter) = sum(dist.^2);
end

withinnss = nan(k,1);
for ic = 1:k
    withinnss(ic) = sum(sum((data(cluster==ic,:)-mean(data(cluster==ic,:),1)).^2,2),1);
end
withinnss = sum(withinnss);

wnull = mean(log(withinnss_null));
w = log(withinnss);
s = std(log(withinnss_null))*sqrt(1+1/niter);
end

%%
function f = fk(data,cluster,kmax)
n = size(data,1);
f = nan(kmax,1);
for ik = 1:kmax
    if ik>1
        Sprev = S;
    end
    S = 0;
    for jc = 1:ik
       S = S+sum(sum((data(cluster(:,ik)==jc,:)-mean(data(cluster(:,ik)==jc,:),1)).^2,2));
    end
    if ik==1
        a = 0;
        f(ik) = 1;
    else
        if ik==2
            a = 1-3/(4*n);
        else
            a = a+(1-a)/6;
        end
        if S==0
            f(ik) = 1;
        else
            f(ik) = S/(a*Sprev);
        end
    end
end
end