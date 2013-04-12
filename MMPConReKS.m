function selectedLeafs = MMPConReKS(data, y, treeStruct, varNames, selectedLeafs, MRoptions, dataInfo)

allNodes=treeStruct.groupMembers.Node;
selectedLeafs.varNames = varNames;

for i=1:length(selectedLeafs.ID)
    i
    % if the ancestors haven't been initialized
    if ~isfield(selectedLeafs,'selectedNodeID') || length(selectedLeafs.selectedNodeID) < i
        % based on the selected leafs, grab all their ancestors
        % arrange their ancestors from youngest to oldest
        selectedLeafs.selectedNodeID{i} = flipud(find(cellfun(@(x)ismember(selectedLeafs.ID(i),x),allNodes)));
        % initialize collapsedData to empty
        selectedLeafs.collapsedData{i} = [];
    end
    
    selectedLeafs.flag{i} = zeros(1,length(selectedLeafs.selectedNodeID{i})); % overall flag
    
    j=1;
    csData = data(:,setdiff(selectedLeafs.ID, selectedLeafs.ID(i))); % need to change to accomodate maxK, may or may not need
    csIndex = 2:(size(csData,2)+1);  % need to change to accomodate maxK, may or may not need
    xData = data(:,selectedLeafs.ID(i));
    collapsed_xIndex = 1;
    
    
    while j <= length(selectedLeafs.selectedNodeID{i}) % keep climbing the tree while we haven't reached the top or haven't been interrupted
        
        if  size(selectedLeafs.collapsedData{i},2)<j % if we haven't created collapsed data yet
            % collapse these genes
            toCollapse = treeStruct.groupMembers.get(selectedLeafs.selectedNodeID{i}(j));
            selectedLeafs.size{i}(j) = length(toCollapse);
            
            % collapse these nodes using collapsing method of choice
            collapsedData = collapseNode(data(:,toCollapse), MRoptions.method, y);
            selectedLeafs.collapsedData{i}(:,j) = collapsedData;
        else % otherwise, use existing collapased data to save time
            collapsedData = selectedLeafs.collapsedData{i}(:,j);
        end
        
        % as long as collapsedData is not invalid (possible with FA)
        if ~isnan(collapsedData)
            % option 1: test (x',T|[cs\x]), calculate parent node significance, using CIT of choice, to
            % modify to include up to 5 subgroups
            [selectedLeafs.xTcsP{i}(j) selectedLeafs.xTcsstat{i}(j) flag uniModelFit selectedLeafs.dev1cs{i}(j) selectedLeafs.dev2xcs{i}(j)]= testIndLogistic(y, [collapsedData csData], collapsed_xIndex, csIndex, dataInfo);
            
            % option 2: test (x,T|x'), calculate parent node significance, using CIT of choice,
            [selectedLeafs.xTxpP{i}(j) selectedLeafs.xTxstat{i}(j) flag uniModelFit selectedLeafs.dev3xp{i}(j) selectedLeafs.dev4xxp{i}(j)]= testIndLogistic(y, [xData collapsedData], 1, 2, dataInfo);
            
            % option 3: test (x,T|[cs\x x']), calculate parent node significance, using CIT of choice,
            %[selectedLeafs.xTcsxpP{i}(j) selectedLeafs.xTcsxpstat{i}(j) flag uniModelFit selectedLeafs.dev5csxp{i}(j) selectedLeafs.dev6csxxp{i}(j)]= testIndLogistic(y, [xData collapsedData csData], 1, [2 csIndex], dataInfo);
            
            % if parent node significance pass both option 1 and option 2 criteria, accept and move up to test older ancester
            % if (selectedLeafs.xTcsP{i}(j) < MRoptions.threshold1 ) & (selectedLeafs.xTxpP{i}(j) > MRoptions.threshold2) % option 1 and 2
            % if selectedLeafs.xTcsP{i}(j) < MRoptions.threshold1 % option 1
            if selectedLeafs.xTxpP{i}(j) > MRoptions.threshold2 % option 2
            % if selectedLeafs.xTcsxpP{i}(j) > MRoptions.thresh % option 3
                                
                selectedLeafs.flag{i}(j)=1; % option 1 and 2
                j=j+1;
                % else stop here and go on to test other leafs % option 1
            else
                break,end
            
        else % if collapsedData is invalid, tests are unnecessary and flag is too
            selectedLeafs.flag{i}(j)=NaN;
            j=j+1;
        end
        
    end
end