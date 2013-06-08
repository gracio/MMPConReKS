function [cv,clusters] = calcStability(cv)
% calculate stability over all cv runs

% first pick out the selected variables
for k=1:length(cv.MR)
    for i=1:length(cv.MR{k}.selectedLeafs)
        for m=1:size(cv.MR{k}.selectedLeafs{i}.groupVarSelected,1)
            for n=1:size(cv.MR{k}.selectedLeafs{i}.groupVarSelected,2)
                for j=1:size(cv.MR{k}.selectedLeafs{i}.groupVarSelected,3)
                    clusters{i,k,m,n}{j} = getFieldByThresh(cv.MR{k}.selectedLeafs{i},cv.MR{k},m,n,j,'toCollapse');
                end
            end
        end
    end
end

% calculate group stability
for i=1:length(cv.MR{k}.selectedLeafs)
    for m=1:size(cv.MR{k}.selectedLeafs{i}.groupVarSelected,1)
        for n=1:size(cv.MR{k}.selectedLeafs{i}.groupVarSelected,2)
            cumVal2 = [];
            for k1=1:length(cv.MR)
                for k2 = k1+1:length(cv.MR)
                    [val1(k1,k2) mi{k1,k2} mj{k1,k2}]=cluster_matching(clusters{i,k1,m,n},clusters{i,k2,m,n}); % maximum bipartite matching
                    val2(k1,k2) = length(mi{k1,k2})/(length(clusters{k1})+length(clusters{k2})-length(mi{k1,k2})); % tanimoto score on top
                    cumVal2 = [cumVal2;val2(k1,k2)];
                end
            end
            %cv.stability{i}(m,n) = mean(val1(val1>0));
            cv.stability{i}(m,n) = mean(cumVal2);
        end
    end
end

% calculate baseline (single) stability
cumStab = [];
for k1=1:length(cv.MR)
    for k2 = k1+1:length(cv.MR)
        stab_single(k1,k2) = length(intersect(cv.MR{k1}.ID,cv.MR{k2}.ID))/(length(union(cv.MR{k1}.ID,cv.MR{k2}.ID)));
        cumStab = [cumStab;stab_single(k1,k2)];
    end
end

% taking average of baseline(single) stability
cv.singleStability = mean(cumStab);