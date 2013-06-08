function plotAccuracyHeatmap(cv)
% plot heatmap of group var selection accuracy across cv runs. calibrate
% the colors according to overall min and max

limLow =1;
limHigh = 0;

for i=1:length(cv.groupMMPCoptions.method)
    limLow = min(limLow, min(min( mean(cv.accuracy{i},3))));
    limHigh = max(limHigh, max(max( mean(cv.accuracy{i},3))));
end

range = limHigh-limLow;
limHigh = min(limHigh + range/1.5, 1);
limLow = max(limLow - range/4,0);

limLow = min(limLow, mean(cv.singleAccuracy));
limHigh = max(limHigh, mean(cv.singleAccuracy));

for i=1:length(cv.groupMMPCoptions.method)
    figure('color',[1 1 1])
    imagesc(mean(cv.accuracy{i},3),[limLow, limHigh])
    colormap('hot')
    axis square
    title(['simulated data: SVM CV accuracy rate, ' cv.groupMMPCoptions.method{i}],'fontsize',14)
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
    
    saveas(gcf,[cv.basicRunPath 'AccuracyHeatMap ' cv.groupMMPCoptions.method{i}],'fig')
end

% plot accuracy for single variables: SVM cv error
figure('color',[1 1 1])
imagesc(mean(cv.singleAccuracy),[limLow, limHigh])
colormap('hot')
xlabel(['single variable baseline accuracy: ' num2str(mean(cv.singleAccuracy))],'fontsize',14)
saveas(gcf,[cv.basicRunPath 'AccuracyHeatMap single'],'fig')