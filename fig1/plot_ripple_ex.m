clearvars; clc; close all;

load('D:\OneDrive\1.allen-andermann\Totalinfo.mat','session_metric');
sessionList = session_metric.session_id(session_metric.immobile_period>100 & session_metric.distance_bw_ml_probes>2000);
nSession = length(sessionList);

rippletype = {'global';'medial';'lateral'};
iS = 1;
load([sdir(sessionList(iS)),'_ripples.mat'],'spontaneous_CA1_ripple','CA1_ripple_classified');
ilist = [1,63,8];
fHandle = figure('PaperUnits','Centimeters','PaperPosition',[2 2 20/8*3 20/8*2]);
[~,sortidx] = sort(cellfun(@(x) x(3),spontaneous_CA1_ripple.ccf_coordinate));
lfp = cellfun(@(x) [x(:,1),nanmean(x(:,2:4),2)],spontaneous_CA1_ripple.lfp(sortidx),'UniformOutput',false);
for iRp = 1:3
   axes('Position',axpt(3,1,iRp,1));
   hold on;
   for iP = 1:length(lfp)
       if iP==1
           clr = 'r';
       elseif iP==length(lfp)
           clr = 'b';
       else
           clr = 'k';
       end   
      plot(lfp{iP}(:,1),lfp{iP}(:,2)+(iP-1)*10^(-3),'Color',clr); 
   end
   
   if isnan(ilist(iRp))
   i = 0;
   satisfied = 0;
   while satisfied==0
       i = i+1;
       xlim([-0.05 0.1]+CA1_ripple_classified.(rippletype{iRp})(i,1));
       satisfied = input('satisfied?');
   end
   ilist(iRp) = i;
   else
      i = ilist(iRp); 
           xlim([-0.05 0.1]+CA1_ripple_classified.(rippletype{iRp})(i,1));
   end
   if iRp==1
   plot([-0.04 0.01]+CA1_ripple_classified.(rippletype{iRp})(i,1),[-0.7 -0.7]*10^(-3),'k');
   text(-0.05+CA1_ripple_classified.(rippletype{iRp})(i,1),-0.4*10^(-3),'0.5s','FontSize',6,'Color','k')
   end
   set(gca,'Box','off','TickDir','out','FontSize',7,'LineWidth',0.35,...
       'XTick',[],'YTick',([1 6]-1)*10^(-3),'YTickLabel',{'M';'I'},'YLim',[-1 6.5]*10^(-3));
   if iRp>1
       set(gca,'YTickLabel',[]);
   end
end
cd('D:\OneDrive - University of California, San Francisco\figures\allen\fig1');
print(fHandle,'-depsc','-painters','example_lfp.ai');