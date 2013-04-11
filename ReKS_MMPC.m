resultPath = '/Volumes/Macintosh HD 2/work/benosWork/projects/MMPConReKS/results/'
dataName = 'BR Tumor ER option 1&2';
if ~isdir([resultPath dataName])
    mkdir([resultPath dataName])
end

git add ReKS_MMPC.m
git commit -m "run check commit version"
[status,cmdout]  = system('git rev-list --max-count=1 HEAD')
mkdir([resultPath dataName '/' num2str(cmdout(1:10)) ])

diary([resultPath dataName '/' num2str(cmdout(1:10)) '/runDiary.log'])

fprintf(['code version ' num2str(cmdout) '\n'])

dataStruct = BR_Tumor;

data = dataStruct.data';
treeStruct = dataStruct.ReKS.treeStruct;
varNames = dataStruct.genes;

y = BR_Tumor.ER;


% MMPConReKS options
MRoptions.threshold1 =  0.05;
MRoptions.threshold2 =  0.05;
method = {'centroid','medoid','CA','PCA','FA'};
MRoptions.method = method{i};
MRoptions.plot = 1;

MRfile

% dataInfo
dataInfo.yType = 'categorical';
dataInfo.yCounts = max(y) + 1;
dataInfo.dataType = 'continuous';
dataInfo.dataCounts = nan;
options.test = 'testIndLogistic';
% MMPC options
options.threshold = 0.05;
options.maxK = 5;


% run single MMPC
if ~exist('selectedLeafs_single')
    [selectedLeafs_single.ID selectedLeafs_single.pvalues selectedLeafs_single.statistics] = MMPC(y, data, dataInfo, options);
end

if ~exist('selectedLeafs')
    for i=1:length(method)
        fprintf(['running ' method{i} '...\n'])
        [selectedLeafs{i}, toplot{i}] = MMPConReKS(data, y,treeStruct, thresh1, thresh2, method{i}, dataName, varNames, selectedLeafs_single, dataInfo);
    end
end

% build classifier for every combination of selected variable/groups
    SVMStruct=[]; SVMindex=[];
    for i=1:length(method)
        % extract all combinations of significant nodes
        varComboN = cellfun(@length,selectedLeafs{i}.parentP)-1;
        
        iszero = find(varComboN==0);
        nonzero = find(varComboN>0);
        if ~isempty(nonzero)
            combos2loopthrough = makePermutationN(varComboN(varComboN>0)); 
        end
        
        SVMindex{i}(:,iszero) = zeros(size(combos2loopthrough,1),length(iszero));
        SVMindex{i}(:,nonzero) = combos2loopthrough;
        
        % loop through combinatinos of significant nodes
        singletonData = [];
        singletonData = data(:,selectedLeafs_single.ID(iszero));
        %size(combos2loopthrough,1)
        for j=1:size(combos2loopthrough,1)
            [i j]
            if ~isempty(nonzero)
            groupData=[];
            for k=1:length(nonzero)
                groupData(:,k) = selectedLeafs{i}.collapsedData{nonzero(k)}(:, combos2loopthrough(j,k));
            end
            svmX_train = [singletonData groupData];
            else
                svmX_train = singletonData;
            end
            svmY_train = y;
            SVMStruct{i}{j} = svmtrain(svmX_train,svmY_train);
        end
    end

save(dataName,'selectedLeafs*', 'toplot', 'SVMStruct', 'SVMindex')

% plotting size against significance
linecolor = 'bgkmc';

for i=1:length(selectedLeafs_single.ID)
    maxx=0;
    figure
    for j=1:length(method)
        xx = [1 selectedLeafs{j}.size{i}(1:length(selectedLeafs{j}.parentP{i}))];
        yy = [selectedLeafs_single.pvalues(selectedLeafs_single.ID(i)) selectedLeafs{j}.parentP{i}];
        yy(yy==0) = 10^-15;
        loglog(xx ,yy ,['.-' linecolor(j)], 'markersize', 15)
        hold on
        if selectedLeafs{j}.size{i}(length(selectedLeafs{j}.parentP{i})) > maxx
            maxx = selectedLeafs{j}.size{i}(length(selectedLeafs{j}.parentP{i}));
        end
    end
    xlabel('node size')
    ylabel('log (IndLogistic) P value')
    title(['node size vs. discriminating significance of ' varNames(selectedLeafs{j}.ID(i))])
    loglog([1 maxx],[thresh thresh],'r-')
    legend([method 'thresh'])
    saveas(gcf,[batchName dataName '-' varNames(selectedLeafs{j}.ID(i))],'fig')
end


%case study on the first selecte nodes:
i=1
% figure out a roughly 2:3 ratio for arranging subplots
figure
for j=1:18
    subplot(3,6,j)
    toCollapse = BR_Tumor.ReKS.treeStruct.groupMembers.get(testSelectedNodes{i}(j));
    collapsedData = collapseNode(data(:,toCollapse), method, y);
    plot(data([find(BR_Tumor.ER==0); find(BR_Tumor.ER==1)],toCollapse),'.')
    hold on
    plot(collapsedData([find(BR_Tumor.ER==0); find(BR_Tumor.ER==1)]),'k.','markersize',15)
    title(['level from leaf : ' num2str(j)])
    ylim([-8 8])
end
figure
for j=19:length(testSelectedNodes{i})
    subplot(3,5,j-18)
    toCollapse = BR_Tumor.ReKS.treeStruct.groupMembers.get(testSelectedNodes{i}(j));
    collapsedData = collapseNode(data(:,toCollapse), method, y);
    plot(data([find(BR_Tumor.ER==0); find(BR_Tumor.ER==1)],toCollapse),'.')
    hold on
    plot(collapsedData([find(BR_Tumor.ER==0); find(BR_Tumor.ER==1)]),'k.','markersize',15)
    title(['level from leaf : ' num2str(j)])
    ylim([-8 8])
end