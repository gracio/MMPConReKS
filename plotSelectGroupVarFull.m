function [cumPval, nodePvalSizes] = plotSelectGroupVarFull(selectedLeafs,treeStruct,title, runPath,varargin)

visible = 1; % default values unless otherwise specified
quick = 0; % default values unless otherwise specified

if nargin ==5
    visible = varargin{1};
elseif nargin == 6
    quick = varargin{2};
end

% get a list of all p-values
cumPval{1} = [];
cumPval{2} = [];

for j=1:length(selectedLeafs.selectedNodeID)
    
    cumPval{1} = [cumPval{1} selectedLeafs.xTcsP{j}];
    cumPval{2} = [cumPval{2} selectedLeafs.xTxpP{j}];
    
    for p=1:length(selectedLeafs.xTcsP{j})
        nodePvalSizes{j}(p,:) = [log10(selectedLeafs.xTcsP{j}(p)) log10(selectedLeafs.xTxpP{j}(p)) selectedLeafs.size{j}(p)];
    end
    
    if quick ~=1
        % plot on the tree
        node2display = [];
        displayContent = [];
        node2display = [node2display; selectedLeafs.selectedNodeID{j}];
        for p=1:length(selectedLeafs.xTcsP{j})
            tp = {['{' num2str(log10(selectedLeafs.xTcsP{j}(p)),2) '/' num2str(log10(selectedLeafs.xTxpP{j}(p)),2) '/' num2str(selectedLeafs.size{j}(p)) '}']};
            displayContent = [displayContent; tp];
        end
        
        h = plotTreeSelectTreeLabel(treeStruct, node2display, displayContent);
        xlabel([title ' variable ' num2str(selectedLeafs.ID(j)) ' single var p-value: ' num2str(selectedLeafs.pvalues(selectedLeafs.ID(j)))],'fontsize',14,'color',[0 0 0])
        saveas(gcf, [runPath title '_FullTree_SelVar_' num2str(j)], 'fig')
        if visible ~= 1
            close(h)
        end
    end
    
end