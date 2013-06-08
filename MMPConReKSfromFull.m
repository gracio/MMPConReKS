function selectedLeafs = MMPConReKSfromFull(selectedLeafs,thresh1,thresh2)

for m=1:length(thresh1)
    for n=1:length(thresh2)
        for j=1:length(selectedLeafs.ID)
            
            selectedLeafs.groupVarThreshGrids(m,n,1) = thresh1(m);
            selectedLeafs.groupVarThreshGrids(m,n,2) =  thresh2(n);
            flag = (selectedLeafs.xTcsP{j} < thresh1(m)) & (selectedLeafs.xTxpP{j} > thresh2(n));
            % find the first occurence of 0
            for i=1:length(flag)
                if flag(i)==0
                    e=i;
                    break,end
            end
            e=e-1;
            selectedLeafs.groupVarSelected(m,n,j) = e;
            
        end
    end
end