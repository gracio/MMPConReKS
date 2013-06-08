function plotMaxClustHeatmap(cv)
% plot heatmap of group var selection max cluster size across cv runs. calibrate
% the colors according to overall min and max

limLow =1;
limHigh = 0;

for i=1:length(cv.groupMMPCoptions.method)
    limHigh = max(limHigh, max(max(cv.maxClustSize{i})));
end

limHigh = limHigh + limHigh/5;

for i=1:length(cv.groupMMPCoptions.method)
    figure('color',[1 1 1])
    imagesc(cv.maxClustSize{i},[limLow, limHigh])
    colormap('hot')
    axis square
    title(['simulated data: SVM CV max cluster size, ' cv.groupMMPCoptions.method{i}],'fontsize',14)
    colorbar
    % set the x axis to the right label
    for xx=1:length(cv.groupMMPCoptions.thresh2)
        labX{xx} = num2str(log10(cv.groupMMPCoptions.thresh2(xx)),2);
    end
    % set the y axis to the right label
    set(gca,'XtickLabel',labX)
    for yy=1:length(cv.groupMMPCoptions.thresh1)
        labY{yy} = num2str(log10(cv.groupMMPCoptions.thresh1(yy)),2);
    end
    set(gca,'YtickLabel',labY)
    xlabel('log10 p-value threshold for ind(x,T|x")','fontsize',14)
    ylabel('log10 p-value threshold for dep(x",T|cs\x)','fontsize',14)
    
    saveas(gcf,[cv.basicRunPath 'maxClustSizeHeatMap ' cv.groupMMPCoptions.method{i}],'fig')
end