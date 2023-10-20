
function [pairid,noisecorr,avepupil,avespeed,nbin] = noisecorr_nostim(spikeTime,win,binsize,unit_id,rippletime,...
    iState,pupiltime,pupilsize,speedtime,speed,speedlimit)

spiketmp = cellfun(@(x) x(x>=win(1) & x<=win(2)),spikeTime,'UniformOutput',false);
spkhist = cell2mat(cellfun(@(x) histcounts(x,win(1):binsize:win(2)),spiketmp,'UniformOutput',false));

rippletmp = rippletime(rippletime>=win(1) & rippletime<=win(2));
ripplehist = histcounts(rippletmp,win(1):binsize:win(2));

[~,~,bin] = histcounts(pupiltime,win(1):binsize:win(2));
pupilhist = cellfun(@(x) nanmean(pupilsize(bin==x)),num2cell(1:floor(diff(win)/binsize)));

[~,~,bin] = histcounts(speedtime,win(1):binsize:win(2));
speedhist = cellfun(@(x) nanmean(speed(bin==x)),num2cell(1:floor(diff(win)/binsize)));

switch iState
    case 1 % every state & no ripple
        incondition = find(ripplehist==0);
    case {2,3} % immobile & no ripple & low-high level of pupil
        pupilmobile = nanmean(pupilhist(abs(speedhist)>speedlimit));
        inconditiontemp = abs(speedhist)<speedlimit & ripplehist==0;
        if iState==2
            incondition = find(inconditiontemp & pupilhist<0.5*pupilmobile);
        elseif iState==3
            incondition = find(inconditiontemp & pupilhist>=0.5*pupilmobile);
        end        
    case 4 % mobile
        incondition = find(abs(speedhist)>speedlimit);
end

if ~isempty(incondition)
    spkhist = spkhist(:,incondition);
    
    nCell = length(spikeTime);
    out = ones(nCell,nCell);
    out = triu(out);
    out = out(:);
    
    x = repmat(unit_id,1,length(unit_id));
    y = repmat(unit_id',length(unit_id),1);
    
    pairid = [x(:),y(:)];
    pairid = pairid(~out,:);
    
    noisecorr = corr(spkhist');
    noisecorr = noisecorr(:);
    noisecorr = noisecorr(~out);
    
    nbin = length(incondition);
    
    avepupil = nanmean(pupilhist(incondition));
    avespeed = nanmean(speedhist(incondition));
else
    nbin = 0;
    [avepupil,avespeed] = deal(NaN);
    noisecorr = NaN(size(pairid,1),1);
end
end