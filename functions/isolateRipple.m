function [isolatedRippletime,isolatedindex] = isolateRipple(rippletimes, mininterval, intervaltype)
%ISOLATERIPPLE: select isolated ripples whose interval to neighboring
%ripples are larger than mingap.
%
%rippletimes: nx2 matrix. Each row corresponds to each ripple incident, and
%first and second column corresponsd to start and end time of each rppple,
%respectively.
%mininterval: minimum interval between isolated ripple events.
%intervaltype: 'onoff','on'

rippleinitial = rippletimes;
rippletimes = sortrows(rippletimes,1);
if strcmp(intervaltype,'on')
    deltarippletimes = [[nan;diff(rippletimes(:,1))],...
        [diff(rippletimes(:,1));nan]];
elseif strcmp(intervaltype,'onoff')
    deltarippletimes = [rippletimes(:,1)-[nan;rippletimes(1:end-1,2)],...
        [rippletimes(2:end,1);nan]-rippletimes(:,2)];
end
isolatedRippletime = rippletimes(sum(deltarippletimes<mininterval,2)==0,:);
isolatedindex = ismember(rippleinitial(:,1),isolatedRippletime(:,1));