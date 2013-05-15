function toplot = plotSelectedLeafs(selectedLeafs, treeStruct, runPath, runName, MRoptions)
% currently plots the bar graph of MMPConReKS

% specify the colors of :(< threshold, > threshold, NaN) in bar graph
colors = {[0 154 205]/255;[139 28 98]/255;[238 238 0]/255};

toplot = [];
for i=1:length(selectedLeafs.ID)
    toplot(i,:) = [sum(selectedLeafs.flag{i}==1) sum(selectedLeafs.flag{i}==0) sum(isnan(selectedLeafs.flag{i}))];
end

hasNaN = logical(sum(toplot(:,3)));
legendLabels = {['ind(x";T|[cs\x] p < ' num2str(MRoptions.threshold1) ' & ind(x;T|x" p > ' num2str(MRoptions.threshold2)],['ind(x";T|[cs\x] p >= ' num2str(MRoptions.threshold1)  ' & ind(x;T|x" p <= ' num2str(MRoptions.threshold2)], 'NaN'}; % option 1 & 2
%legendLabels = {['ind(x";T|[cs\x] p < ' num2str(MRoptions.threshold1)],['ind(x";T|[cs\x] p >= ' num2str(MRoptions.threshold1)], 'NaN'}; % option 1
%legendLabels = {['ind(x;T|x" p > ' num2str(MRoptions.threshold2)],['ind(x;T|x" p <= ' num2str(MRoptions.threshold2)], 'NaN'}; % option 2

figure


hBar = bar([toplot; 0 0 0],'stacked')
set(hBar,{'FaceColor'},colors);
legend(legendLabels(sum(toplot,1)>0));


title(runName)
xlabel('selected variables')
ylabel('level counting from leaf')
set(gca,'XTickLabel',selectedLeafs.varNames(selectedLeafs.ID)')
saveas(gcf,[runPath runName],'fig')

%%%
node2display = [];
contents = [];
for i=1:length(selectedLeafs.selectedNodeID)
    node2display = [node2display; selectedLeafs.selectedNodeID{i}(selectedLeafs.flag{i}>0)];
    contents = [contents;  repmat(selectedLeafs.varNames(selectedLeafs.ID(i)),sum(selectedLeafs.flag{i}>0), 1)];
end


% if there are multiplle contents to be displayed at the same node,
% concatenate them

uqNode = unique(node2display);
k=1;
for i=1:length(uqNode)
    if sum(ismember(node2display,uqNode(i)))>1
        dup{k} = find(ismember(node2display,uqNode(i)));
        k=k+1;
    end
    
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
if ~isempty(node2display)
plotTreeSelectTreeLabel(treeStruct.nodeID, node2display, contents)
end
saveas(gcf,[runPath 'clusterVarSelectedVarOnTree'],'fig')
