function collapsedData = collapseNode(data2collapse,method, y)
% Summarize the group of variables using method of choice
% methods are: FA, PCA, CA, centroid, medoid, along second dimension

if size(data2collapse,2)>1

switch method
    case 'centroid'
        collapsedData = mean(data2collapse,2);
    case 'medoid'
        [m,medoidI]=min(sum(squareform(pdist(data2collapse'))));
        collapsedData = data2collapse(:,medoidI);
    case 'FA' % significantly slow
        try
        [loadings,specVar,t,stats,collapsedData] = factoran(data2collapse,1);
        catch
            collapsedData = NaN;
        end
    case 'PCA'
        [coefs,collapsedData,var,tw]=princomp(data2collapse);
        collapsedData = collapsedData(:,1);
    case 'CA'
        [A,B,r,U,V,stats] = canoncorr(data2collapse(~isnan(y),:),y(~isnan(y)));
        [m,CAI] = max(abs(A));
        collapsedData = data2collapse(:,CAI);
end

else
    collapsedData = data2collapse;
end