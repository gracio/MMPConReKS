function [parents] = getParents(treeStruct, nodeID)

parents = [];
nodeIDD = nodeID;

while nodeIDD~=1
    nodeIDD = treeStruct.nodeID.getparent(nodeIDD);
    parents = [parents nodeIDD];
end
