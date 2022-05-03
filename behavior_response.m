clc; close all; clearvars;
load('D:\heejeong\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
load('D:\heejeong\OneDrive\1.allen-andermann\tag.mat');

sessionList = session_metric.session_id(session_metric.immobile_period>100 &...
    session_metric.global_ripple_number>0);
nS = length(sessionList);
for iS = 12:nS
    iS
    try
       load([sdir(sessionList(iS)),'_stimulus_evoked.mat'],'fr_behavior');
    end
    if exist('fr_behavior','var')
        clear fr_behavior
        continue;
    end    
    
    [selec_run,pval_selec_run,selec_pupil,pval_selec_pupil,...
        r_corr_run,pval_corr_run,r_corr_pupil,pval_corr_pupil] = beh_metric(sessionList(iS));
    
    fr_behavior.run.selec.idx = selec_run;
    fr_behavior.run.selec.pval = pval_selec_run;
    fr_behavior.run.corr.r = r_corr_run;
    fr_behavior.run.corr.pval = pval_corr_run;
    
    fr_behavior.pupil.selec.idx = selec_pupil;
    fr_behavior.pupil.selec.pval = pval_selec_pupil;
    fr_behavior.pupil.corr.r = r_corr_pupil;
    fr_behavior.pupil.corr.pval = pval_corr_pupil;
    
    save([sdir(sessionList(iS)),'_stimulus_evoked.mat'],'fr_behavior','-append');
end

function [selec_run,pval_selec_run,selec_pupil,pval_selec_pupil,...
    r_corr_run,pval_corr_run,r_corr_pupil,pval_corr_pupil] = beh_metric(sessionid)

load([sdir(sessionid),'_ripples.mat'],'spontaneous_win','spontaneous_anal_win');
load([sdir(sessionid),'_Events.mat'],'running_speed','pupil_data');
load([sdir(sessionid),'_cellTable.mat'],'T');

%% running
[spkhist,spkhist_r] = deal({});
for iW = 1:size(spontaneous_anal_win,1)
    if diff(spontaneous_anal_win(iW,:))>1
   spkhist(:,iW) = cellfun(@(x) histcounts(x,spontaneous_anal_win(iW,1):1:spontaneous_anal_win(iW,2)),...
       T.spike_time,'UniformOutput',false);
    end
   if iW==1
       if spontaneous_anal_win(iW,1)>spontaneous_win(1)
           spkhist_r = [spkhist_r,cellfun(@(x) histcounts(x,spontaneous_win(1):1:spontaneous_anal_win(iW,1)),...
               T.spike_time,'UniformOutput',false)];
       end
   elseif iW==size(spontaneous_anal_win,1)
       if spontaneous_anal_win(iW,2)<spontaneous_win(2)
           spkhist_r = [spkhist_r,cellfun(@(x) histcounts(x,spontaneous_anal_win(iW,2):1:spontaneous_win(2)),...
               T.spike_time,'UniformOutput',false)];
       end
   else
       spkhist_r = [spkhist_r,cellfun(@(x) histcounts(x,spontaneous_anal_win(iW-1,2):1:spontaneous_anal_win(iW,1)),...
           T.spike_time,'UniformOutput',false)];
   end
end
spkhist = mat2cell(cell2mat(spkhist),ones(size(spkhist,1),1),sum(cellfun(@length,spkhist(1,:))),1);
spkhist_r = mat2cell(cell2mat(spkhist_r),ones(size(spkhist_r,1),1),sum(cellfun(@length,spkhist_r(1,:))),1);
selec_run = [cellfun(@mean,spkhist_r)-cellfun(@mean,spkhist)]./sqrt(cellfun(@std,spkhist_r).^2+cellfun(@std,spkhist).^2);
[~,pval_selec_run] = cellfun(@(x,y) ttest2(x,y),spkhist,spkhist_r);

run_velocity = running_speed.velocity;
run_time = running_speed.time;
inrun = run_time>=spontaneous_win(1) & run_time<=spontaneous_win(2);
run_velocity = abs(run_velocity(inrun));
run_time = run_time(inrun);

spkhist_run = cellfun(@(x) histc(x,run_time),T.spike_time,'UniformOutput',false);
spkhist_run = cellfun(@(x) x*(1/nanmean(diff(run_time))),spkhist_run,'UniformOutput',false);

[rtmp,ptmp] = cellfun(@(x) corrcoef(run_velocity,x,'Rows','complete'),spkhist_run,'UniformOutput',false);
r_corr_run = cellfun(@(x) x(1,2),rtmp);
pval_corr_run = cellfun(@(x) x(1,2),ptmp);

%% pupil
pupil_size = (pupil_data.pupil_height/2).*(pupil_data.pupil_width/2)*pi;
inpupil = pupil_data.time>=spontaneous_win(1) & pupil_data.time<=spontaneous_win(2);
pupil_size = pupil_size(inpupil);
pupil_time = pupil_data.time(inpupil);

[~,sortIdx] = sort(pupil_size);
low25 = sortIdx(1:round(0.25*length(sortIdx)));
high25 = sortIdx(length(sortIdx)-round(0.25*length(sortIdx)):end);

spkhist_pupil = cellfun(@(x) histc(x,pupil_time),T.spike_time,'UniformOutput',false);
spkhist_pupil = cellfun(@(x) x*(1/nanmean(diff(pupil_time))),spkhist_pupil,'UniformOutput',false);

spkhigh = cellfun(@(x) x(high25),spkhist_pupil,'UniformOutput',false);
spklow = cellfun(@(x) x(low25),spkhist_pupil,'UniformOutput',false);
selec_pupil = [cellfun(@mean,spkhigh)-cellfun(@mean,spklow)]./sqrt(cellfun(@std,spkhigh).^2+cellfun(@std,spklow).^2);
[~,pval_selec_pupil] = cellfun(@(x,y) ttest2(x,y),spkhigh,spklow);

[rtmp,ptmp] = cellfun(@(x) corrcoef(pupil_size,x,'Rows','complete'),spkhist_pupil,'UniformOutput',false);
r_corr_pupil = cellfun(@(x) x(1,2),rtmp);
pval_corr_pupil = cellfun(@(x) x(1,2),ptmp);

end

