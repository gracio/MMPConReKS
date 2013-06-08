function cvData = myCV(k, label, repeats, varargin)
% perform balanced cv partition, making sure al label kinds are equally
% distributed, with designated number of repeats
% label is double integers

if ~isempty(varargin)
    cvData = varargin{1};
end

% get all unique labels
uqLabel = unique(label);
uqLabel = uqLabel(~isnan(uqLabel)); % remove all NaN's from consideration

% make cv partitions
for n=1:repeats % for each repeat
    shuffled_id =[];
    % make random permuations of all unique labels for this repeat
    for i=1:length(uqLabel)
        t = find(label == uqLabel(i));
        shuffled_id{i} = t(randperm(length(t)));
        shuffled_id{i} = reshape(shuffled_id{i},1,length(shuffled_id{i}));
    end
    
    for m=1:k % for each fold
        if m < k % first m-1 CV partitions
            
            tp = [];
            for i=1:length(uqLabel)
                tp = [tp shuffled_id{i}( (floor( length(shuffled_id{i})/k )*(m-1)+1) : (floor( length(shuffled_id{i})/k )*m ))];
            end
            
            cvData.testID{n,m} = tp;
            cvData.trainID{n,m} = setdiff(find(~isnan(label)),cvData.testID{n,m}); % train data is the rest, but don't include NaNs
            
            cvData.testSize(n,m) = length(cvData.testID{n,m});
            cvData.trainSize(n,m) = length(cvData.trainID{n,m});
            
        else % special treatment for the last CV partition
            
            tp = [];
            for i=1:length(uqLabel)
                tp = [tp shuffled_id{i}( floor(length(shuffled_id{i})/k)*(m-1)+1 : length(shuffled_id{i}) )];
            end
            
            cvData.testID{n,m} = tp;
            cvData.trainID{n,m} = setdiff(find(~isnan(label)),cvData.testID{n,m}); % train data is the rest, but don't include NaNs
            
            cvData.testSize(n,m) = length(cvData.testID{n,m});
            cvData.trainSize(n,m) = length(cvData.trainID{n,m});
            
        end
    end
end

cvData.testID = reshape(cvData.testID,1,k*repeats);
cvData.trainID = reshape(cvData.trainID,1,k*repeats);
cvData.testSize = reshape(cvData.testSize,1,k*repeats);
cvData.trainSize = reshape(cvData.trainSize,1,k*repeats);