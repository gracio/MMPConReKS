function [dataStruct] = MMPConReKS_main(dataStruct,runPath,minClustSize,singleMMPCoptions,groupMMPCoptions,dataInfo,runName)


visible = 0; % turn off display on screen
quick = 1; % do not plot ReKS trees and selected single variables up the tree
lite = 1; % do not generate extra information for ReKS tree: A, S, U, discComp, centroid, dist from top and bottom

% generate ReKS tree
dataStruct = ReKS_main(dataStruct,runPath,minClustSize,dataStruct.varDim,visible,quick,lite); % turn off display on screen to avoid run out of mem

% make sure Y is a column vector
dataStruct.Y = reshape(dataStruct.Y,length(dataStruct.Y),1);

tic
% run single MMPC
[selectedLeafs_single.ID selectedLeafs_single.pvalues selectedLeafs_single.statistics] = MMPC(dataStruct.Y, dataStruct.data, dataInfo, singleMMPCoptions);
%plotSingleDiagnostics(data,y,varNames,treeStruct,selectedLeafs_single,runPath);
fprintf(['single variable selection finished...\n'])
toc

% put the selected varaibles in the data structure
dataStruct.ID = selectedLeafs_single.ID;

tic
for i=1:length(groupMMPCoptions.method)
    
    MRoptions.method = groupMMPCoptions.method{i};
    fprintf(['running ' MRoptions.method '...\n'])
    
    % first, run the full tree to level specified. if no level is
    % specified, run all the way to the top
    fprintf('running full tree to level \n')
    MRoptions.thresh1 =  1.1; % most lenient
    MRoptions.thresh2 =  -0.1; % most lenient
    MRoptions.levelLimit = groupMMPCoptions.levelLimit;
    dataStruct.selectedLeafs{i} = MMPConReKS(dataStruct.data,dataStruct.Y,dataStruct.ReKS.treeStruct,dataStruct.genes, selectedLeafs_single, MRoptions, dataInfo);
    
    % get the p-value distribution and plot the tree if quick is off
    [cumPval, dataStruct.selectedLeafs{i}.nodePvalSizes] = plotSelectGroupVarFull(dataStruct.selectedLeafs{i},dataStruct.ReKS.treeStruct.numDescendants,[runName ' ' groupMMPCoptions.method{i}],runPath,visible,quick); % turn off display on screen to avoid run out of mem
    
    % identify p-value range automatically
    fprintf(['plotting p-value distribution for ' MRoptions.method '...\n'])
    figure('color',[1 1 1])
    for c=1:length(cumPval)
        d=log10(cumPval{c});
        d(d==-Inf)=-15;
        subplot(length(cumPval),1,c)
        hist(d,40)
        PvalMin{c} = min(d(d~=-15))
        PvalMax{c} = max(d)
        maxer = ylim;
        maxer = maxer(2);
        hold on
        if c==1
            for m=1:length(groupMMPCoptions.thresh1)
                line(repmat(log10(groupMMPCoptions.thresh1(m)),2,1),[0 maxer])
            end
            title([groupMMPCoptions.method{i} ': distribution of dep(x",T|cs\x) pvalues on the full tree'])
        else
            for m=1:length(groupMMPCoptions.thresh2)
                line(repmat(log10(groupMMPCoptions.thresh2(m)),2,1),[0 maxer])
            end
            title([groupMMPCoptions.method{i} ': distribution of ind(x,T|x") pvalues on the full tree'])
        end
        saveas(gcf,[runPath groupMMPCoptions.method{i} 'p-val distribution'],'fig')
    end
    fprintf(['done plotting p-value distribution for ' MRoptions.method '...\n'])
    
    if groupMMPCoptions.run2top ==1
        % get the results from each thershold
        dataStruct.selectedLeafs{i} = MMPConReKSfromFull(dataStruct.selectedLeafs{i},groupMMPCoptions.thresh1,groupMMPCoptions.thresh2);
        fprintf(['done thresholding for ' MRoptions.method '...\n'])
        
    else
        
        % if user doesn't run the single variables all the way to the top, do
        % as much as needed
        fprintf('running tree for as much as needed \n')
        % obtain the trees at various levels
        MRoptions.thresh1 =  groupMMPCoptions.thresh1; % most lenient
        MRoptions.thresh2 =   groupMMPCoptions.thresh2; % most lenient
        MRoptions.levelLimit = Inf;
        
        dataStruct.selectedLeafs{i} = MMPConReKS(dataStruct.data,dataStruct.Y,dataStruct.ReKS.treeStruct,dataStruct.genes, dataStruct.selectedLeafs{i}, MRoptions, dataInfo);
        %plotSelectedLeafs(selectedLeafs{i,m,n}, treeStruct, runPath, [runName ' ' method{i}], MRoptions);
    end
    
end

dataStruct.groupMMPCoptions = groupMMPCoptions;
dataStruct.singleMMPCoptions = singleMMPCoptions;

dataStruct = rmfield(dataStruct, 'ReKS');
dataStruct = rmfield(dataStruct, 'genes');
dataStruct = rmfield(dataStruct, 'headers');
dataStruct = rmfield(dataStruct, 'sparseSymA');
fprintf('thresholding finished \n')
toc
