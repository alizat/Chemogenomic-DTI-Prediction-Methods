function CV_indices=getCVindices(seed,cv_setting,N,Y)
    % randomly generate indices for the N folds
    CV_indices = ones(numel(Y), 1);
    rng(seed)  % seed for reproducibility of results
    [num_drugs,num_targets] = size(Y);
    if strcmp(cv_setting,'S1')
        len = numel(Y);
    elseif strcmp(cv_setting,'S2')
        len = num_drugs;
    elseif strcmp(cv_setting,'S3')
        len = num_targets;
    end
    rng('default')
    rng(seed);
    rand_ind = randperm(len);
    for i=1:N
        % leave out random drug-target pairs
        if strcmp(cv_setting,'S1')
            test_ind = rand_ind((floor((i-1)*len/N)+1:floor(i*len/N))');
            left_out = test_ind;

        % leave out random entire drugs
        elseif strcmp(cv_setting,'S2')
            left_out_drugs = rand_ind((floor((i-1)*len/N)+1:floor(i*len/N))');
            test_ind = zeros(length(left_out_drugs),num_targets);
            for j=1:length(left_out_drugs)
                curr_left_out_drug = left_out_drugs(j);
                test_ind(j,:) = ((0:(num_targets-1)) .* num_drugs) + curr_left_out_drug;
            end
            test_ind = reshape(test_ind,numel(test_ind),1);
            left_out = left_out_drugs;

        % leave out random entire targets
        elseif strcmp(cv_setting,'S3')
            left_out_targets = rand_ind((floor((i-1)*len/N)+1:floor(i*len/N))');
            test_ind = zeros(num_drugs,length(left_out_targets));
            for j=1:length(left_out_targets)
                curr_left_out_target = left_out_targets(j);
                test_ind(:,j) = (1:num_drugs)' + ((curr_left_out_target-1)*num_drugs);
            end
            test_ind = reshape(test_ind,numel(test_ind),1);
            left_out = left_out_targets;
        end
        left_out = left_out(:);
        test_ind = test_ind(:);

        CV_indices(test_ind) = i;
    end
end