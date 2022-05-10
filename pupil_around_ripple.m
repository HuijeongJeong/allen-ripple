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

win = [-20 20];
binsize = 0.1;
bin = win(1):binsize:win(2);
rippletype = {'medial';'lateral';'global'};
pupil_ripple = cell(nS,3);

for iS = 1:nS
    iS
    load([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_anal_win',...
        'spontaneous_win','CA1_ripple_classified');
    load([sdir(sessionList(iS)),'_Events.mat'],'pupil_data','running_speed');
    
    pupil_size = pupil_data.pupil_height/2.*pupil_data.pupil_width/2*pi;
    running_win = [spontaneous_anal_win(1:end-1,2), spontaneous_anal_win(2:end,1)];
    if spontaneous_anal_win(1,1)>spontaneous_win(1)
        running_win = [spontaneous_win(1) spontaneous_anal_win(1); running_win];
    end
    if spontaneous_anal_win(end,2)<spontaneous_win(2)
        running_win = [running_win; spontaneous_anal_win(end,2) spontaneous_win(2)];
    end
    
    inrunning = false(length(pupil_size),1);
    for iW = 1:size(running_win,1)
        inrunning = inrunning | [pupil_data.time>=running_win(iW,1) & pupil_data.time<=running_win(iW,2)];
    end
    pupil_running = nanmean(pupil_size(inrunning));
    pupil_size_norm = pupil_size/pupil_running;
    
    for iRp =1:3
        rippletime = CA1_ripple_classified.(rippletype{iRp})(:,1);
        intime = cellfun(@(x) pupil_data.time>=x+win(1)-binsize & pupil_data.time<=x+win(2)+binsize,...
            num2cell(rippletime),'UniformOutput',false);
        pupil_ripple{iS,iRp} = cell2mat(cellfun(@(x,y) interp1(pupil_data.time(x)-y,pupil_size_norm(x),bin),...
            intime,num2cell(rippletime),'UniformOutput',false));
    end
end

m = cellfun(@nanmean,pupil_ripple,'UniformOutput',false);
ct = {'r';'b';'k'};

fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 6 2.5]);
for i = 1:2
    axes('Position',axpt(2,10,i,1:8,[],[0.2 0.05]));
hold on;
for iRp = 1:3
    mm = nanmean(cell2mat(m(:,iRp)));
    ms = nanstd(cell2mat(m(:,iRp)))/sqrt(nS);
    fill([bin, flip(bin)],[mm+ms flip(mm-ms)],ct{iRp},'EdgeColor','none');
    plot(bin,mm,'color',ct{iRp});
end
plot([0 0],[0 1],'k:');
alpha(0.1);
xlabel('Time from ripple (s)')
ylabel('Normalized pupil');
set(gca,'Box','off','TickDir','out','FontSize',5,'LineWidth',0.35);
if i==1
   ylim([0 1]); 
else
    ylim([0.35 0.65]);
end
end
cd('D:\OneDrive - University of California, San Francisco\figures\allen');
print(fHandle,'-dtiff','-r600','pupil_around_ripple.tif');

iS
function sessionDir = sdir(sessionid)

sessionDir = ['D:\OneDrive\1.allen-andermann\session_',num2str(sessionid),...
    '\session_',num2str(sessionid)];
end
