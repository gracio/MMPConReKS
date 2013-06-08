function out = getFieldByThresh(selectedLeafs,dataStruct,mm,nn,jj,fieldName)
% mm is a vector(or single) of threshold1's to loop over. nn is a vector(or single)
% of thershold2's to loop over. jj is the index of the variable we wish to
% extract information for

for m_it=1:length(mm) % loop over all members of mm's
    m=mm(m_it);
    
    for n_it=1:length(nn) % loop over all members of nn's
        n=nn(n_it);
        
        for j_it=1:length(jj) % loop over all members of jj's
            j=jj(j_it);
            
            % get the nodeID of mth threh1, nth thresh2, and ith selected variable
            id = selectedLeafs.groupVarSelected(m,n,j);
            
            if id==0
                switch fieldName
                    case 'selectedNodeID'
                        out(m_it,n_it,j_it) = selectedLeafs.selectedNodeID{j}(1); % output the lowest node if one wants single var
                    case 'collapsedData'
                        if selectedLeafs.varDim == 1
                            out{m_it,n_it,j_it} = dataStruct.data(selectedLeafs.ID(j),:)';
                        else
                            out{m_it,n_it,j_it} = dataStruct.data(:,selectedLeafs.ID(j));
                        end
                    case 'size'
                        out(m_it,n_it,j_it) = 1;
                    case 'xTcsP'
                        out(m_it,n_it,j_it) = NaN;
                    case 'xTxpP'
                        out(m_it,n_it,j_it) = NaN;
                    case 'toCollapse'
                        out{m_it,n_it,j_it} = selectedLeafs.ID(j);
                end
                
            else
                switch fieldName
                    case 'selectedNodeID'
                        out(m_it,n_it,j_it) = selectedLeafs.selectedNodeID{j}(id);
                    case 'collapsedData'
                        out{m_it,n_it,j_it} = selectedLeafs.collapsedData{j}(:,id);
                    case 'size'
                        out(m_it,n_it,j_it) = selectedLeafs.size{j}(id);
                    case 'xTcsP'
                        out(m_it,n_it,j_it) = selectedLeafs.xTcsP{j}(id);
                    case 'xTxpP'
                        out(m_it,n_it,j_it) = selectedLeafs.xTxpP{j}(id);
                    case 'toCollapse'
                        out{m_it,n_it,j_it} = selectedLeafs.toCollapse{j}{id};
                end
            end
        end
    end
end

if iscell(out) & length(out)==1
    out = out{1};
end
