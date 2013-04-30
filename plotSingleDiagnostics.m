function plotSingleDiagnostics(data,y,varNames,treeStruct,selectedLeafs_single,savePath)
% take a look at all the selected single variables


% plot clustergram
if length(selectedLeafs_single.ID)>1
t=selectedLeafs_single.pvalues(selectedLeafs_single.ID);
for i=1:length(selectedLeafs_single.ID)
    rowlab{i} = [varNames{selectedLeafs_single.ID(i)} ' (' num2str(selectedLeafs_single.pvalues(selectedLeafs_single.ID(i)), '%0.1e') ')'];
end

c = clustergram(data([find(y==0) ;find(y==1)],selectedLeafs_single.ID)','rowlabels',rowlab,'columnlabels',[zeros(1,sum(y==0)) ones(1,sum(y==1))])
c.plot
saveas(gcf,[savePath 'singleVarClustergram'],'fig')

order = c.RowLabels;
for i=1:length(order)
    t = order{i};
    startDelete = strfind(t,' ');
    t(startDelete:length(t))=[];
neworder(i) = strmatch(t,varNames(selectedLeafs_single.ID));
end
neworder = fliplr(neworder);
else
    neworder = 1:length(selectedLeafs_single.ID);
end

% plot correlation
figure
imagesc(corr(data([find(y==0) ;find(y==1)],selectedLeafs_single.ID(neworder))))
saveas(gcf,[savePath 'singleVarCorrelation'],'fig')

% plot the labels
figure
imagesc([zeros(1,sum(y==0)) ones(1,sum(y==1))])
saveas(gcf,[savePath 'singleVarYlabel'],'fig')

allNodes=treeStruct.groupMembers.Node;
% find out where the selected single variables are on the tree
for i=1:length(selectedLeafs_single.ID)
    t = flipud(find(cellfun(@(x)ismember(selectedLeafs_single.ID(i),x),allNodes)));
    node2display(i) = t(1);
end

for i=1:length(node2display)
    contents(i) = varNames(selectedLeafs_single.ID(i));
end

% if there are multiplle contents to be displayed at the same node,
% concatenate them

uqNode = unique(node2display);
k=1;
for i=1:length(uqNode)
    if sum(ismember(node2display,uqNode(i)))>1
        dup{k} = find(ismember(node2display,uqNode(i)));
    end
    k=k+1;
end

if exist('dup')
    toremove = [];
    k = length(node2display)+1;
    for i=1:length(dup)
        newcontents = contents{dup{i}(1)};
        for j=2:length(dup{i})
            newcontents = [newcontents '/' contents{dup{i}(j)}];
        end
        node2display(k) = node2display(dup{i}(1));
        contents{k} = newcontents;
        toremove = [toremove dup{i}];
        k=k+1;
    end
    node2display(toremove) = [];
    contents(toremove) = [];
end

plotTreeSelectTreeLabel(treeStruct.nodeID, node2display, contents)
saveas(gcf,[savePath 'singleVarSelectedVarOnTree'],'fig')
