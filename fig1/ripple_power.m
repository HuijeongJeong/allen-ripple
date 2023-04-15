clearvars; clc; close all;

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

frequency = [100 250];
bandOrder = 4;

win = [-0.5 0.5];
bin = 0.001;
rippletype = {'medial','lateral','global_m','global_l'};

alignednss = cell(nS,6,4);
probe_ccf = cell(nS,1);
for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_Events.mat'],'probe');
    load([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_CA1_ripple',...
        'CA1_ripple_classified','spontaneous_win');
    
    Fs = probe.lfp_sampling_rate;
    [~,probeidx] = sort(cellfun(@(x) x(2),spontaneous_CA1_ripple.ccf_coordinate));
    probe_ccf{iS} = cell2mat(spontaneous_CA1_ripple.ccf_coordinate(probeidx));
    for iP = 1:size(spontaneous_CA1_ripple,1)
        data = spontaneous_CA1_ripple.lfp{probeidx(iP)};
        Fs_probe = round(Fs(probe.id==spontaneous_CA1_ripple.probe_id(probeidx(iP))));
        if sum(isnan(data(:,2)))>0
            nanidx = isnan(data(:,2));
            nanstart = find(diff(nanidx)==1);
            nanend = find(diff(nanidx)==-1);
            analwin = [[1;nanend(:)+1],[nanstart(:);size(data,1)]];
            nss_temp = NaN(size(data,1),size(data,2)-1);
            flfp = NaN(size(data,1),size(data,2));
            for in = 1:size(analwin,1)
                flfp(analwin(in,1):analwin(in,2),:) =...
                    bz_Filter(data(analwin(in,1):analwin(in,2),:),...
                    'passband',frequency,'order',bandOrder,'nyquist',Fs_probe/2); %filter lfp
                nss_temp(analwin(in,1):analwin(in,2),:) = NormalizedSquaredSignal_HJ(flfp(...
                    analwin(in,1):analwin(in,2),:));
            end
        else
            flfp = bz_Filter(data,'passband',frequency,'order',bandOrder,...
                'nyquist',Fs_probe/2); %filter lfp
            nss_temp = NormalizedSquaredSignal_HJ(flfp);
        end
        nss = mean(nss_temp,2);
        lfp_time = data(:,1);
        for iRp = 1:4
            [time,alignednss{iS,iP,iRp}] = alignlfp2event(nss,lfp_time,...
                CA1_ripple_classified.(rippletype{iRp})(:,1),[-1 1],Fs(iP),bin);
        end
    end
end
%%
ct = cbrewer('div','RdBu',7);
ct(3:5,:) = [];
endidx = cellfun(@(x) size(x,1),probe_ccf);
titleList = {'dCA1 ripple';'iCA1 ripple';'Global ripple'};

close all;
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20*0.45 20*0.14]);
for iRp = 1:3
    axes('Position',axpt(3,1,iRp,1,axpt(1,10,1,1:9)));
    hold on;
    for i = [1,4]
        switch i
            case {1,2}
                j = repmat(i,nS,1);
            case 3
                j = endidx-1;
            case 4
                j = endidx;
        end
        data = cell(nS,1);
        for iS = 1:nS
            if iRp<3
                data{iS} = nanmean(alignednss{iS,j(iS),iRp},2)';
            else
                if rem(i,2)==1
                    data{iS} = nanmean(alignednss{iS,j(iS),1},2)'; 
                    % aligned to onset time of ripples detected at dorsal ca1
                else
                    data{iS} = nanmean(alignednss{iS,j(iS),2},2)';
                    % aligned to onset time of ripples detected at intermeidate ca1
                end
            end
        end
        m = mean(cell2mat(data));
        s = std(cell2mat(data))/sqrt(size(data,1));
        fill([time flip(time)],[m+s flip(m-s)],ct(i,:),'EdgeColor','none');
        plot(time,m,'Color',ct(i,:));
    end
    plot([0 0],[-1 12],'k:');
    set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,...
        'XTick',-0.3:0.3:0.3,'YTick',0:4:12);
    title(titleList{iRp});
    if iRp==1
        ylabel({'Norm. ripple'; 'power (z-score)'});
    else
        set(gca,'YTickLabel',[]);
    end
    xlabel('Time from ripple onset (s)');
    xlim([-0.3 0.3]);
    ylim([-1 12]);
end
cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig1')
print(fHandle,'-depsc','-painters','rippl_power_probes.ai');
%%
distance_to_d = cellfun(@(x) sqrt(sum((x-repmat(x(1,:),size(x,1),1)).^2,2)),...
    probe_ccf,'UniformOutput',false);
distance_to_i = cellfun(@(x) sqrt(sum((x-repmat(x(end,:),size(x,1),1)).^2,2)),...
    probe_ccf,'UniformOutput',false);
normdist = cellfun(@(x,y) x./(x+y),distance_to_d,distance_to_i,'UniformOutput',false);
normdist = cellfun(@(x) [x; NaN(6-size(x,1),1)],normdist,'UniformOutput',false);

[peaksize,peaktime] = deal(NaN(nS,4,6));
for iS = 1:nS
    for iRp = 1:4
        for i = 1:6
            if isempty(alignednss{iS,i,iRp})
                continue;
            end
            [peaksize(iS,iRp,i),peaktime(iS,iRp,i)] = max(nanmean(alignednss{iS,i,iRp}(intime,:),2));
        end
    end
end

close all
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20/3 20*0.14]);
for iRp = 1:2
    axes('Position',axpt(2,1,iRp,1,axpt(10,10,2:10,1:9)));
    hold on;
    x = cell2mat(normdist')';
    y = squeeze(peaksize(:,iRp,:));
    beta = glmfit(x(:),y(:));
    [r,p] = corr(x(:),y(:),'Rows','Complete');
    scatter(x(:),y(:),2,[0.6 0.6 0.6],'filled')
    scatter(x(x==0),y(x==0),2,ct(1,:),'filled')
    scatter(x(x==1),y(x==1),2,ct(end,:),'filled');
    plot([0 1],[0 1]*beta(2)+beta(1),'k');
    ylim([0 16]);
    xlim([-0.1 1.1]);
    set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,...
        'XTick',0:0.5:1,'YTick',0:5:15);
    if iRp==1
        ylabel({'Peak norm.'; 'ripple power (z-score)'});
        xlabel('Norm. distance from dCA1 probe');
    else
        set(gca,'YTickLabel',[]);
    end
    title(titleList{iRp});
end
print(fHandle,'-depsc','-painters','ripple_power_distance.ai');

%%
function [time,alignedlfp] = alignlfp2event(lfp,lfptime,eventtime,win,Fs,bin)
Fs = round(Fs);
[~,idx] = cellfun(@(x) min(abs(lfptime-x)),num2cell(eventtime));
alignedlfp = cell2mat(cellfun(@(x) lfp([round(win(1)*Fs):round(win(2)*Fs)]+x),...
    num2cell(idx),'UniformOutput',false)');
if size(alignedlfp,1)~=round(diff(win)*Fs+1)
    alignedlfp = alignedlfp';
end
alignedlfp = interp1([win(1):1/Fs:win(2)]',alignedlfp,[win(1):bin:win(2)]);
time = win(1):bin:win(2);
end

