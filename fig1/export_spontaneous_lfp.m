clc; clearvars; close all;

load('D:\heejeong\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\heejeong\OneDrive\1.allen-andermann\tag.mat');
load('D:\heejeong\OneDrive\1.allen-andermann\ripple_glm\ripple_glm_cluster.mat');

bandOrder = 4;

[~,idx] = ismember(unit_id.vis,tag.info.unit_id);
session_id = tag.info.session_id(idx);

unitid = unit_id.vis;
sessionList = unique(session_id);
nS = length(sessionList);

 for iS = 2:nS
     iS
     load([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_win','spontaneous_CA1_ripple');  
     load([sdir(sessionList(iS)),'_Events.mat'],'probe');
     
     probe_id = spontaneous_CA1_ripple.probe_id;
     
     for iP = 1:length(probe_id)
         load(['D:\heejeong\Documents\AllenSDK_data\session_',num2str(sessionList(iS)),...
             '\probe_',num2str(probe_id(iP)),'_lfp.mat'],'lfp_3','lfp_4','channel_id','time');
         
         intime = cellfun(@(x) x>=spontaneous_win(1)-10 & x<=spontaneous_win(2)+10,...
             time(3:4),'UniformOutput',false);
         inchannel = ismember(channel_id,spontaneous_CA1_ripple.channel_id{iP});
         
         lfp = cell2mat(cellfun(@(x,y) x(inchannel,y),{lfp_3,lfp_4},intime,'UniformOutput',false));
         time = cell2mat(cellfun(@(x,y) x(y),time(3:4),intime,'UniformOutput',false));
         
         spontaneous_CA1_ripple.lfp{iP} = [time',lfp'];
     end
     save([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_CA1_ripple','-append');
 end
 
 
function dir = sdir(sessionnum)
dir = ['D:\heejeong\OneDrive\1.allen-andermann\session_',num2str(sessionnum),'\session_',num2str(sessionnum)];
end
