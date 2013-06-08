% acquire (real) data

load /Volumes/Macintosh' HD 2'/work/benosWork/projects/MMPConReKS/data/test_real/BR_Tumor.mat
cv.runName =  'BR Tumor ER CV'
cv.data = BR_Tumor.data
cv.Y = BR_Tumor.ER
cv.genes = BR_Tumor.genes;
cv.headers = BR_Tumor.headers
cv.varDim = 1

% load BNresults_sameBN.mat
% cv.data = trainData{1}.data
% cv.Y = trainData{1}.Y
% cv.genes = trainData{1}.genes;
% cv.headers = trainData{1}.headers
% cv.varDim = 1
% cv.runName = 'test simulate'

clearvars -except cv

resultPath = '/Volumes/Macintosh HD 2/work/benosWork/projects/MMPConReKS/results/'; % folder where I store results

cv.k = 5; % k-fold CV
cv.repeats = 2; % CV repeats
cv.minClustSize = 5; % minimum cluster size for use in ReKS
cv.thresh2plot = []; % don't specify which thresholds to plot - default of diagonal is used

cv.groupMMPCoptions.method = {'centroid','medoid','CA','PCA'};
%groupMMPCoptions.method = {'centroid','medoid','CA','PCA','FA'};
cv.groupMMPCoptions.thresh1 = [10^-9 10^-8 10^-7 10^-6 10^-5 10^-4 5*(10^-4) (10^-3) 5*(10^-3) (10^-2) ];
cv.groupMMPCoptions.thresh2 = [10^-4 5*(10^-4) 10^-3 10^-2 5*(10^-2) 10^-1 2*(10^-1) 4*(10^-1) 6*(10^-1) 8*(10^-1)];
cv.groupMMPCoptions.run2top = 0; % either run to the top of the tree, then do thresholding, or do each threshold at a time
cv.groupMMPCoptions.levelLimit = 5; % either run to the top of the tree, then do thresholding, or do each threshold at a time

cv.singleMMPCoptions.threshold = 0.05;
cv.singleMMPCoptions.maxK = 5;
cv.singleMMPCoptions.test = 'testIndLogistic';

% dataInfo
cv.dataInfo.yType = 'categorical';
cv.dataInfo.yCounts = max(cv.Y) + 1;
cv.dataInfo.dataType = 'continuous';
cv.dataInfo.dataCounts = nan;
cv.dataInfo.varDim = cv.varDim; %perform variable selection along first dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate Cross Validation datasets
cv = myCV(cv.k, cv.Y, cv.repeats, cv);
fprintf('generated cross validation data sets\n')

% generate the commitID/dateTime folder for this result
[cv.basicRunPath, cv.commitID] = resultVCfolderHelper(cv.runName, resultPath);

% for each cross validation run, do:
for k=1:length(cv.trainID)
    
    % generate folder for this run
    runPath = [cv.basicRunPath 'cv' num2str(k) '/'];
    mkdir(runPath)
    
    % log the parameters and run and commit reftag
    diary([runPath 'runDiary.log'])
    fprintf(['code version ' num2str(cv.commitID) '\n'])
    tic
    
    % generate training data for this run
    if cv.varDim == 1
        trainCV.data = cv.data(:,cv.trainID{k});
    else
        trainCV.data = cv.data(cv.trainID{k},:);
    end
    trainCV.Y = cv.Y(cv.trainID{k});
    trainCV.genes = cv.genes;
    trainCV.headers = cv.genes(cv.trainID{k});
    trainCV.varDim = cv.varDim;
    fprintf(['cv run ' num2str(k) ' size of training data: ' num2str(size(trainCV.data)) '\n'])
    
    % generate test data for this run
    if cv.varDim == 1
        testCV.data = cv.data(:,cv.testID{k});
    else
        testCV.data = cv.data(cv.testID{k},:);
    end
    testCV.Y = cv.Y(cv.testID{k});
    testCV.varDim = cv.varDim;
    fprintf(['cv run ' num2str(k) ' size of testing data: ' num2str(size(testCV.data)) '\n'])
    
    % run main file to perform variable selections across thresholds
    [cv.MR{k}] = MMPConReKS_main(trainCV,runPath,cv.minClustSize,cv.singleMMPCoptions,cv.groupMMPCoptions,cv.dataInfo,cv.runName);
    
    % train SVM
    cv.MR{k} = trainSVMmodel(cv.MR{k});
    
    % test SVM and get accuracy
    cv.testResult{k} = testSVMmodel(cv.MR{k},testCV);
    
    % create a table of single variables selected over runs
    cv.singleTable(k,1:length(cv.MR{k}.ID)) = cv.MR{k}.ID;
    
    % create a table of group variables selected over runs.....only works
    % for simulated data. pending.
    
    toc
    % end log files
    diary
    
end

% calculate stability
fprintf(['calculating stability \n'])
cv = calcStability(cv);

fprintf(['creating summary statistics \n'])
m=length(cv.groupMMPCoptions.thresh1);
n=length(cv.groupMMPCoptions.thresh2);
% create summary statistics averaged over CV runs:
for k=1:length(cv.MR)
    % summarize single variable selection cv accuracy
    cv.singleAccuracy(k) = cv.testResult{k}.singleAccuracy;
    
    for i=1:length(cv.groupMMPCoptions.method)
        % summarize group variable selection cv accuracy
        cv.accuracy{i}(:,:,k) = cv.testResult{k}.accuracy{i};
        % summarize group variable mean and max cluster sizes
        cv.meanClustSize{i} = mean(getFieldByThresh(cv.MR{k}.selectedLeafs{i},cv.MR{k},1:m,1:n,1:length(cv.MR{k}.ID),'size'),3);
        cv.maxClustSize{i} = max(getFieldByThresh(cv.MR{k}.selectedLeafs{i},cv.MR{k},1:m,1:n,1:length(cv.MR{k}.ID),'size'),[],3);
    end
end

% save this mat file
clear i k m n ans
save([cv.basicRunPath 'cv.mat'], '-v7.3')

% generate plots for all CV runs:
fprintf(['plotting results \n'])
cv.thresh2plot = plotYYYAccStabClust(cv); % yyy line plots. using default option: plotting diagonal
plotAccuracyHeatmap(cv) % heatmap: accuracy + single
plotStabilityHeatmap(cv) % heatmap: stability  + single
plotMeanClustHeatmap(cv) % heatmap: mean size
plotMaxClustHeatmap(cv) % heatmap: max size
