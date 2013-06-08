function thresh2plot = plotYYYAccStabClust(cv,varargin)
% generate YYY line plot

thresh1 = cv.groupMMPCoptions.thresh1;
thresh2 = cv.groupMMPCoptions.thresh2;
nThresh = min(length(thresh1),length(thresh2));

% figure out which thresholds to plot
if ~isempty(varargin)
    thresh2plot = varargin{1}; % if user supplies which thresholds to plot
else
    % use defalt: plot across the diagonal
    thresh2plot(:,1) = 1:nThresh; % do the min to accomdate cases where thresh 1 and thresh 2 are different sizes
    thresh2plot(:,2) = length(thresh2)-(1:nThresh)+1;
end

% generate threshold labels (log10)
for j=1:nThresh
    lab{j} = [num2str(log10(thresh1(thresh2plot(j,1))),1) '/' num2str(log10(thresh2(thresh2plot(j,2))),2)];
end

% extract the statistics based on the sampling: row is sample, column is
% method i
for i=1:length(cv.groupMMPCoptions.method)
    for j=1:nThresh
        mm = thresh2plot(j,1);
        nn = thresh2plot(j,2);
        
        % extract accuracy
        tp = mean(cv.accuracy{i},3);
        acc(j,i)= tp(mm,nn);
        % extract stability
        stab(j,i)= cv.stability{i}(mm,nn);
        % extract max cluster size
        clustMax(j,i)= cv.maxClustSize{i}(mm,nn);
        clustMean(j,i)= cv.meanClustSize{i}(mm,nn);
    end
end

% now start plotting
for i=1:length(cv.groupMMPCoptions.method)
  
    % plot 1st y axis: accuracy and mean accuracy(dotted), 2nd y axis:
    % stability and mean stability(dotted, 3rd: max cluster size
    [ax,hlines,figh] = plotyyy(1:nThresh,[acc(:,i) repmat(mean(cv.singleAccuracy),nThresh,1)],1:nThresh,[stab(:,i) repmat(cv.singleStability,nThresh,1)],1:nThresh,clustMax(:,i),{'accuracy','stability','max cluster size'});
    
    title(cv.groupMMPCoptions.method{i},'fontsize',14,'fontWeight','bold')
    xlabel('log10 p-value thresholds','fontsize',14,'fontWeight','bold')
    set(figh,'color',[1 1 1])
    
    set(get(ax(1),'Ylabel'),'String','Accuracy','fontsize',14,'fontWeight','bold')
    set(get(ax(2),'Ylabel'),'String','Stability','fontsize',14,'fontWeight','bold')
    set(get(ax(3),'Ylabel'),'String','Max Cluster Size','fontsize',14,'fontWeight','bold')

    set(ax(1),'xticklabel',[])
    set(ax(2),'xticklabel',[])
    set(ax(1),'xticklabel',lab,'fontsize',10)
    
    % set accuracy y axis limits and ticks
    trueMin = min(min(acc));
    trueMax = max(max(acc));
    [betterMin, betterMax, betterInc] = getBetterLimits(min(trueMin,mean(cv.singleAccuracy)),max(trueMax,mean(cv.singleAccuracy)),10,[0 1]);
    set(ax(1),'ylim',[betterMin betterMax])
    set(ax(1),'Box','off')
    set(ax(1),'YTick',[betterMin:betterInc:betterMax])
    
    % set stabliity y axis limits and ticks
    trueMin = min(min(stab));
    trueMax = max(max(stab));
    [betterMin, betterMax, betterInc] = getBetterLimits(min(cv.singleStability,trueMin),max(trueMax,cv.singleStability),10,[0 1]);
    set(ax(2),'ylim',[betterMin betterMax])
    set(ax(2),'Box','off')
    set(ax(2),'YTick',[betterMin:betterInc:betterMax])
    
    % set cluster size y axis limits and ticks
    trueMin = min(min(clustMax));
    trueMax = max(max(clustMax));
    [betterMin, betterMax, betterInc] = getBetterLimits(trueMin,trueMax,10,[0 1000]);
    set(ax(3),'ylim',[betterMin betterMax])
    set(ax(3),'Box','off')
    set(ax(3),'YTick',[betterMin:betterInc:betterMax])
    
    % get rid of extra lines 
    limx1=get(ax(1),'xlim');
    limx3=[limx1(1)   limx1(1) + 1.2*(limx1(2)-limx1(1))];
    limy3=get(ax(3),'YLim');
    cfig = get(gcf,'color');
    line([limx1(2) limx3(2)],[limy3(1) limy3(1)],...
        'Color',cfig,'Parent',ax(3),'Clipping','off');
    axes(ax(2))
    
    % fix position of x label
    xlabh = get(gca,'XLabel');
    set(xlabh,'Position',get(xlabh,'Position') - [0 .01 0])
    
    saveas(gcf, [cv.basicRunPath 'YYYplot ' cv.groupMMPCoptions.method{i} ],'fig')
    %close(figh)
end
