function [selectedLeafs, toplot] = MMPConReKS(data, y, treeStruct, thresh1, thresh2, method, dataName, varNames, selectedLeafs, dataInfo)
% specify the colors of :(< threshold, > threshold, NaN) in bar graph
colors = {[0 154 205]/255;[139 28 98]/255;[238 238 0]/255};

allNodes=treeStruct.groupMembers.Node;

for i=1:length(selectedLeafs.ID)
    i
    % based on the selected leafs, grab all their ancestors
    % arrange their ancestors from youngest to oldest
    if ~isfield(selectedLeafs,'selectedNodeID') || length(selectedLeafs.selectedNodeID) < i
        selectedLeafs.selectedNodeID{i} = flipud(find(cellfun(@(x)ismember(selectedLeafs.ID(i),x),allNodes)));
    end
    
    selectedLeafs.flag{i} = zeros(1,length(selectedLeafs.selectedNodeID{i})); % option 1 & 2
    selectedLeafs.xpTcsflag{i} = zeros(1,length(selectedLeafs.selectedNodeID{i})); % option 1
    selectedLeafs.xTxpflag{i} = zeros(1,length(selectedLeafs.selectedNodeID{i})); % option 2
    %selectedLeafs.xTcsxpflag{i} = zeros(1,length(selectedLeafs.selectedNodeID{i})); % option 3
    
    j=1;
    % only needed for option 1 (and 3?)
    csData = data(:,setdiff(selectedLeafs.ID, selectedLeafs.ID(i))); % option 1 & 3
    csIndex = 2:(size(csData,2)+1);% option 1 & 3
    collapsed_xIndex = 1;% option 1
    xData = data(:,selectedLeafs.ID(i)); % option 2 & 3
    
    while j <= length(selectedLeafs.selectedNodeID{i})
        
        if ~isfield(selectedLeafs,'collapsedData') || (length(selectedLeafs.collapsedData) < i) || size(selectedLeafs.collapsedData{i},2)<j
            % collapse these genes
            toCollapse = treeStruct.groupMembers.get(selectedLeafs.selectedNodeID{i}(j));
            selectedLeafs.size{i}(j) = length(toCollapse);
            [i j]
            % collapse these nodes using collapsing method of choice
            collapsedData = collapseNode(data(:,toCollapse), method, y);
            selectedLeafs.collapsedData{i}(:,j) = collapsedData;
        else
            collapsedData = selectedLeafs.collapsedData{i}(:,j);
        end
        
        if ~isnan(collapsedData)
            % option 1: test (x',T|[cs\x]), calculate parent node significance, using CIT of choice, to
            % modify to include up to 5 subgroups
            [selectedLeafs.parentP{i}(j) selectedLeafs.stat{i}(j) flag uniModelFit selectedLeafs.dev1cs{i}(j) selectedLeafs.dev2xcs{i}(j)]= testIndLogistic(y, [collapsedData csData], collapsed_xIndex, csIndex, dataInfo);
            
            % option 2: test (x,T|x'), calculate parent node significance, using CIT of choice,
            [selectedLeafs.xTxpP{i}(j) selectedLeafs.xTxstat{i}(j) flag uniModelFit selectedLeafs.dev3xp{i}(j) selectedLeafs.dev4xxp{i}(j)]= testIndLogistic(y, [xData collapsedData], 1, 2, dataInfo);
            
            % option 3: test (x,T|[cs\x x']), calculate parent node significance, using CIT of choice,
            %[selectedLeafs.xTcsxpP{i}(j) selectedLeafs.xTcsxpstat{i}(j) flag uniModelFit selectedLeafs.dev5csxp{i}(j) selectedLeafs.dev6csxxp{i}(j)]= testIndLogistic(y, [xData collapsedData csData], 1, [2 csIndex], dataInfo);
            
            % for marking purposes
            if (selectedLeafs.parentP{i}(j) < thresh1 )% option 1
                selectedLeafs.xpTcsflag{i}(j) = 1;% option 1
            end
            if (selectedLeafs.xTxpP{i}(j) > thresh2) % option 2
                selectedLeafs.xTxpflag{i}(j)=1; % option 2
            end
            
            % if parent node significance pass both option 1 and option 2 criteria, accept and move up to test older ancester
            if (selectedLeafs.parentP{i}(j) < thresh1 ) & (selectedLeafs.xTxpP{i}(j) > thresh2) % option 1 and 2                
                selectedLeafs.flag{i}(j)=1; % option 1 and 2
                j=j+1;
                % else stop here and go on to test other leafs % option 1
            else
                break,end
            
        else
            selectedLeafs.flag{i}(j)=NaN;
            j=j+1;
        end
        
        %                     if selectedLeafs.xTxpP{i}(j) > thresh % option 2
        %                 selectedLeafs.xTxpflag{i}(j)=1; % option 2
        %                 j=j+1; % option 2
        %                 % else stop here and go on to test other leafs
        %             else % option 2
        %                 break,end % option 2
        %
        %         else % option 2
        %             selectedLeafs.xTxpflag{i}(j)=NaN; % option 2
        %             j=j+1; % option 2
        %         end % option 2
        
        
        %             if selectedLeafs.xTcsxpP{i}(j) > thresh % option
        %                 selectedLeafs.xTcsxpflag{i}(j)=1; % option 3
        %                 j=j+1; % option 3
        %                 % else stop here and go on to test other leafs
        %             else % option 3
        %                 break,end % option 3
        %
        %         else % option 3
        %             selectedLeafs.xTcsxpflag{i}(j)=NaN; % option 3
        %             j=j+1; % option 3
        %         end % option 3
    end
end

toplot = [];
for i=1:length(selectedLeafs.ID)
    toplot(i,:) = [sum(selectedLeafs.flag{i}==1) sum(selectedLeafs.flag{i}==0) sum(isnan(selectedLeafs.flag{i}))]; % option 1
    %toplot(i,:) = [sum(selectedLeafs.xTxpflag{i}==1) sum(selectedLeafs.xTxpflag{i}==0) sum(isnan(selectedLeafs.xTxpflag{i}))]; % option 2
    %toplot(i,:) = [sum(selectedLeafs.xTcsxpflag{i}==1) sum(selectedLeafs.xTcsxpflag{i}==0) sum(isnan(selectedLeafs.xTcsxpflag{i}))]; % option 3
end

hasNaN = logical(sum(toplot(:,3)));

figure
if hasNaN
    hBar = bar(toplot,'stacked')
    set(hBar,{'FaceColor'},colors);
    legend(['p < ' num2str(thresh)],['p >= ' num2str(thresh)], 'NaN')
else
    toplot(:,3)=[];
    hBar = bar(toplot,'stacked')
    set(hBar,{'FaceColor'},colors(1:2));
    legend(['ind(x";T|[cs\x] p < ' num2str(thresh1) ' & ind(x;T|x" p > ' num2str(thresh2)],['ind(x";T|[cs\x] p >= ' num2str(thresh1)  ' & ind(x;T|x" p <= ' num2str(thresh2)])
end

title([dataName ' -' method])
xlabel('selected variables')
ylabel('level counting from leaf')
set(gca,'XTickLabel',varNames(selectedLeafs.ID)')
saveas(gcf,[dataName ' -' method],'fig')