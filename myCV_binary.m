function cv_id_toDelete = myCV_binary(cvpartition_m, label)

% make random permuations of both labels
t=find(label==1);
shuffled_id_1 = t(randperm(length(t)));
t=find(label==0);
shuffled_id_0 = t(randperm(length(t)));

for m=1:cvpartition_m
    m

    if m < cvpartition_m % first m-1 CV partitions
        cv_id_toDelete{m} = shuffled_id_1( floor(length(shuffled_id_1)/cvpartition_m)*(m-1)+1 : floor(length(shuffled_id_1)/cvpartition_m)*m ); % pick out the 1's
        cv_id_toDelete{m} = [cv_id_toDelete{m}; shuffled_id_0( floor(length(shuffled_id_0)/cvpartition_m)*(m-1)+1 : floor(length(shuffled_id_0)/cvpartition_m)*m )]; % pick out the 0's

    else % special treatment for the last CV partition
        cv_id_toDelete{m} = shuffled_id_1( floor(length(shuffled_id_1)/cvpartition_m)*(m-1)+1 : length(shuffled_id_1) ); % pick out the 1's
        cv_id_toDelete{m} = [ cv_id_toDelete{m}; shuffled_id_0( floor(length(shuffled_id_0)/cvpartition_m)*(m-1)+1 : length(shuffled_id_0) ) ]; % pick out the 0's

    end
end