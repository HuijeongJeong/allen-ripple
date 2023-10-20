clearvars;clc;close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\OneDrive\1.allen-andermann\tag.mat');
load('D:\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.distance_bw_ml_probes>2000 & session_metric.global_ripple_number>0);
nS = length(sessionList);

fwin = [0.15, 1];

threshold = [0.3, 0.5];
win = [-3, 3];
binsize = 0.01;
binsize_pupil = 0.0333;

psession = nan(nS,5,811);
psession_pupil = nan(nS,2700);
[phaseripple,magnituderipple] = deal(cell(nS,6));
numvalley = zeros(nS,6);
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_cellTable.mat'],'T');
    load([sdir(sessionList(iS)),'_Events.mat'],'pupil_data');
    load([sdir(sessionList(iS)),'_ripples.mat'],'CA1_ripple_classified','spontaneous_anal_win');
    
    pupil_size = pupil_data.pupil_height/2.*pupil_data.pupil_width/2*pi;
    
    %%
    targets = tag.info.unit_id(tag.area.vis&tag.celltype.rs | tag.area.thalamus);
    in = ismember(T.unit_id,targets);
    
    unitid = T.unit_id(in);
    spiketime = T.spike_time(in);
    
    cidx = zeros(sum(in),1);
    [in,idx] = ismember(unitid,unit_id.vis);
    cidx(in) = cluster_idx.vis{2}(idx(in));    
    
    [in,idx] = ismember(unitid,tag.info.unit_id(tag.area.thalamus));
    cidx(in) = 4;
    %%
    spkhist = cell2mat(cellfun(@(x) histcounts(x,...
        spontaneous_anal_win(1,1)-10:binsize:spontaneous_anal_win(end,2)+10),...
        spiketime,'UniformOutput',false));
    avespkhist = nan(5,size(spkhist,2));
    for iclst = 1:5
        avespkhist(iclst,:) = mean(spkhist(cidx==iclst-1,:));
    end
    time = spontaneous_anal_win(1,1)+binsize/2-10:binsize:spontaneous_anal_win(end,2)-binsize/2+10;
    
    %%
    ptotal = cell(6,1);
    for iwin = 1:size(spontaneous_anal_win,1)
        if diff(spontaneous_anal_win(iwin,:))<30
            continue;
        end
        inwin = time>=spontaneous_anal_win(iwin,1) & time<=spontaneous_anal_win(iwin,2);
        rippleinwin = CA1_ripple_classified.global(CA1_ripple_classified.global(:,1)>=...
            spontaneous_anal_win(iwin,1) & CA1_ripple_classified.global(:,1)<=spontaneous_anal_win(iwin,2),1);
        for iclst = 1:6
            if iclst<6
                if sum(cidx==iclst-1)<5
                    continue;
                end
                [p,f] = pspectrum(avespkhist(iclst,inwin),1/binsize);
                filtered = bandpass(zscore(avespkhist(iclst,inwin)),fwin,1/binsize);
                timeanal = time(inwin);
                f_clst = f;
            else
                inwinanal = pupil_data.time>=spontaneous_anal_win(iwin,1) &...
                    pupil_data.time<=spontaneous_anal_win(iwin,2);
                timeanal = pupil_data.time(inwinanal);
                data = pupil_size(inwinanal);
                if sum(isnan(data))>0
                    data = interp1(timeanal(~isnan(data)),data(~isnan(data)),timeanal);
                    timeanal = timeanal(~isnan(data));
                    data = data(~isnan(data));
                end
                [p,f] = pspectrum(data,1/binsize_pupil);
                filtered = bandpass(zscore(data),fwin,1/binsize_pupil);
                f_pupil = f;
            end
            ptotal{iclst} = [ptotal{iclst};p(f>0.1 & f<=10)'];
            x = hilbert(filtered);
            magnitude = abs(x);
            phase = angle(x);
            
            %%
            invalley = phase<=-pi/2 | phase>=pi/2;
            startvalley = find(diff([nan,invalley(:)'])==1);
            endvalley = find(diff([invalley(:)',nan])==-1);
            
            numvalley(iS,iclst) = numvalley(iS,iclst)+min([length(startvalley),length(endvalley)]);
            %%
            phaseripple{iS,iclst} = [phaseripple{iS,iclst};...
                cellfun(@(x) phase(find(timeanal<x,1,'last')),num2cell(rippleinwin))];
            magnituderipple{iS,iclst} = [magnituderipple{iS,iclst};...
                cellfun(@(x) magnitude(find(timeanal<x,1,'last')),num2cell(rippleinwin))];
        end
    end
    %%
    psession(iS,~cellfun(@isempty,ptotal(1:5)),:) =...
        cell2mat(cellfun(@(x) mean(x,1)/sum(mean(x,1)),ptotal(~cellfun(@isempty,ptotal(1:5))),'UniformOutput',false));
    psession_pupil(iS,:) = mean(ptotal{6},1)/sum(mean(ptotal{6},1));
end
%% eaxmple
clstname = {'Nomod','iAct','dAct','Inh','Thal','Pupil'};

close all
resolution = 2;
avespkconv = conv2(avespkhist,fspecial('Gaussian',[1 5*resolution],resolution),'same');
clr = [[0.6 0.6 0.6];cbrewer('qual','Dark2',4)];
clrrp = {'k','r','b'};
rippletype = {'global','medial','lateral'};

window = [5400, 5500];
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 12 6]);
j = 1;
for i = [4,5,2,3,1]
    inwin = time>=window(1) & time<=window(2);
    filtered = bandpass(zscore(avespkhist(i,inwin)),fwin,1/binsize);
    axes('Position',axpt(1,5,1,j));
   plot(time(inwin),zscore(avespkconv(i,inwin)),'Color',clr(i,:)); 
   hold on;
   plot(time(inwin),filtered,'Color',[0.2 0.2 0.2])
   set(gca,'Box','off','TickDir','out','FontSize',8,'LineWidth',0.35,...
       'XTick',window(1):20:window(2));
   if j==5
      xlabel('Time (s)');
      ylabel({clstname{i};'Norm. firing (z-score)'});
   else
      set(gca,'XTickLabel',[]); 
      ylabel({clstname{i};''});
      if j==1
          for irp = 1:3
             scatter(CA1_ripple_classified.(rippletype{irp})(:,1),...
                 ones(size(CA1_ripple_classified.(rippletype{irp}),1),1)*4,...
                 3,clrrp{irp},'filled');
          end
      end
   end
   xlim(window)
   
   j = j+1;
end
% print(fHandle,'-dtiff','-r600','D:\OneDrive - UCSF\data\allen\revision\example_oscillation_psth.tif');
% print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\data\allen\revision\example_oscillation_psth.ai');
%%
out = sum(isnan(psession(:,:,1)),2)>0;
close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 12 4]);
for i = 1:5
   subplot(1,5,i)
        f = f_clst;
        data = squeeze(psession(~out,i,:));
    imagesc(f(f>0.1 & f<=10),1:nS,data);
    xlim([0.1 2])
    set(gca,'CLim',[0 0.012])
    title(clstname{i});
end
% print(fHandle,'-dtiff','-r600','D:\OneDrive - UCSF\data\allen\revision\bandpower_heatmap.tif');

%%
data_ref = squeeze(psession(:,1,:));
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 12 3]);
for i = 1:4
%     subplot(1,4,i)
    axes('Position',axpt(4,10,i,1:9));
    f = f_clst;
    data = squeeze(psession(:,i+1,:));
    hold on;
    plot(f(f>0.1 & f<=10),data_ref,'Color',[0.8 0.8 0.8])
    plot(f(f>0.1 & f<=10),data,'Color',[1 0.8 0.8])
    plot(f(f>0.1 & f<=10),nanmean(data,1),'r')
    plot(f(f>0.1 & f<=10),nanmean(data_ref,1),'k')    
    xlim([0 2])
    ylim([0 0.04])
    title(clstname{i+1});
    set(gca,'Box','off','TickDir','out','YTick',0:0.01:0.04,'FontSize',8);
    if i==1
        ylabel('Normalized band power');
        xlabel('Frequency (Hz)');
    else
        set(gca,'YTickLabel',[]);
    end
end
% print(fHandle,'-dtiff','-r600','D:\OneDrive - UCSF\data\allen\revision\bandpower.tif');
% print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\data\allen\revision\bandpower.ai');
%%
infrequency = f(f>0.1 & f<=10);
infrequency = infrequency>fwin(1) & infrequency<=fwin(2);
data_ref = mean(squeeze(psession(:,1,infrequency)),2);
%%
for i = 1:4
    %%
    data = mean(squeeze(psession(:,i+1,infrequency)),2);
    [~,p(i)] = ttest(data_ref(~out),data(~out));
end

%%
out = sum(cellfun(@isempty,phaseripple),2)>0;
[histcounts,histbin] = cellfun(@(x) hist(x,25),phaseripple,'UniformOutput',false);
x = cellfun(@(x) x/sum(x),histcounts,'UniformOutput',false);
histbin = histbin{1,1};
clstlist = [4,5,2];
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 9 3]);
for iclst = 1:3
    axes('Position',axpt(3,10,iclst,1:9));
    hold on;
    plot(histbin,cell2mat(x(~out,clstlist(iclst)))','color',[0.8 0.8 0.8]);
    plot(histbin,nanmean(cell2mat(x(~out,clstlist(iclst))),1),'k');
    ylim([0 0.2])
    title(clstname{clstlist(iclst)});
    set(gca,'Box','off','TickDir','out','YTick',0:0.1:0.2,'FontSize',8,...
        'XTick',-pi:pi/2:pi,'XLim',[-pi pi],'XTickLabel',{'-\pi','-1/2\pi',...
        '0','1/2\pi','\pi'},'YTickLabel',0:10:20);
    if iclst==1
        ylabel('% ripple');
        xlabel('Theta');
    end
end
% print(fHandle,'-dtiff','-r600','D:\OneDrive - UCSF\data\allen\revision\phase_at_ripple.tif');
% print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\data\allen\revision\phase_at_ripple.ai');


%%
close all
fraction = cellfun(@(x) sum(x<=-pi/2 | x>=pi/2)/length(x),phaseripple);
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.2 3.5]);
plot(1:3,fraction(~out,clstlist),'Color',[0.8 0.8 0.8]);
hold on;
errorbar(1:3,mean(fraction(~out,clstlist)),std(fraction(~out,clstlist))/sqrt(sum(~out)),'k','CapSize',3)
set(gca,'XTick',1:3,'XTickLabel',clstname(clstlist),'XTickLabelRotation',45,'XLim',[0.5 3.5],'Box','off','TickDir','out');
ylabel('% ripple in valley'); 
% print(fHandle,'-dtiff','-r600','D:\OneDrive - UCSF\data\allen\revision\ripple_in_valley.tif');
print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\data\allen\revision\ripple_in_valley.ai');


%%
fraction = cellfun(@(x) sum(x<=-pi/2 | x>=pi/2),phaseripple)./numvalley;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.2 3.5]);
plot(1:3,fraction(~out,clstlist),'Color',[0.8 0.8 0.8]);
hold on;
errorbar(1:3,mean(fraction(~out,clstlist)),std(fraction(~out,clstlist))/sqrt(sum(~out)),'k','CapSize',3)
set(gca,'XTick',1:3,'XTickLabel',clstname(clstlist),'XTickLabelRotation',45,...
    'XLim',[0.5 3.5],'Box','off','TickDir','out');
ylabel('% valley in ripple'); 
% print(fHandle,'-dtiff','-r600','D:\OneDrive - UCSF\data\allen\revision\valley_in_ripple.tif');
print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\data\allen\revision\valley_in_ripple.ai');

%%
close all
meanphase = cellfun(@circ_mean,phaseripple(~out,:));
data = circ_dist(meanphase(:,4),meanphase(:,5));
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 2.2 3.5]);
hold on;
bar(0.5,mean(data),0.8,'FaceColor',[0.6 0.6 0.6]);
scatter(rand(sum(~out),1)*0.6+0.2,data,3,'k','filled');
errorbar(0.5,mean(data),std(data)/sqrt(sum(~out)),'k','CapSize',3);
ylabel('\Deltaphase(Inh-Thal)');
set(gca,'Box','off','TickDir','out','XTick',[],'FontSize',8,'YLim',[-1/4*pi 1/2*pi],...
    'YTick',[-1/4*pi:1/4*pi:1/2*pi],'YTickLabel',{'-\pi/4','0','\pi/4','\pi/2'},'XLim',[-0.2 1.2]);
print(fHandle,'-depsc','-painters','D:\OneDrive - UCSF\data\allen\revision\delta_phase_inhthal.ai');

[~,p,~,stat] = ttest(data)
