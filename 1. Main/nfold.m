function [auc_res,aupr_res]=nfold(Y,seed,k)
%
% This is a helper function of crossValidation.m. Depending on the
% specified CV setting (or scenario) and supplied "seed", it divides the
% interaction matrix into "nr_fold" folds, performs a cross validation
% experiment and then reports the results (AUPR/AUC) back to
% crossValidation.m.
%
% INPUT:
%  Y:           interaction matrix
%  seed:        seed used for random sampling
%  k:           current repetion of CV experiment
%
% OUTPUT:
%  auc_res:     AUC result
%  aupr_res:    AUPR result
%

    global predictionMethod cv_setting n ds gridSearchMode

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

    if ~gridSearchMode
        fprintf('n-fold experiment start:  \t%s\n',datestr(now));
    end

    AUCs  = zeros(1,n);
    AUPRs = zeros(1,n);
    Y_final = zeros(size(Y));
    % loop over the n folds
    for i=1:n
        % S1: leave out random drug-target pairs
        if strcmp(cv_setting,'S1')
            test_ind = rand_ind((floor((i-1)*len/n)+1:floor(i*len/n))');
            left_out = test_ind;

        % S2: leave out random entire drugs
        elseif strcmp(cv_setting,'S2')
            left_out_drugs = rand_ind((floor((i-1)*len/n)+1:floor(i*len/n))');
            test_ind = zeros(length(left_out_drugs),num_targets);
            for j=1:length(left_out_drugs)
                curr_left_out_drug = left_out_drugs(j);
                test_ind(j,:) = ((0:(num_targets-1)) .* num_drugs) + curr_left_out_drug;
            end
            test_ind = reshape(test_ind,numel(test_ind),1);
            left_out = left_out_drugs;

        % S3: leave out random entire targets
        elseif strcmp(cv_setting,'S3')
            left_out_targets = rand_ind((floor((i-1)*len/n)+1:floor(i*len/n))');
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


        % predict with test set being left out
        y2 = Y;
        y2(test_ind) = 0;   % test set = ZERO
        fprintf('*');       % progress indicator
        y3 = alg_template(y2,predictionMethod,test_ind,left_out); % predict! (y3 is the predicted matrix)


        % compute evaluation metrics based on obtained prediction scores
        AUCs(i)  = calculate_auc (y3(test_ind),Y(test_ind));
        AUPRs(i) = calculate_aupr(y3(test_ind),Y(test_ind));
        if ~gridSearchMode
            fprintf('%.3g\t\t\t\tTIME:    %s\n', AUPRs(i), datestr(now));
        end
        diary off;  diary on;

        Y_final(test_ind) = y3(test_ind);
    end
    
    
    % average the AUCs and AUPRs of the different folds
    auc_res = mean(AUCs);
    aupr_res = mean(AUPRs);


    % print results of the n-fold experiment
    if ~gridSearchMode
        fprintf('n-fold experiment end:  \t%s\n',datestr(now));

        fprintf('\n');
        fprintf('      AUC: %g\n',   auc_res);
        fprintf('     AUPR: %g\n',   aupr_res);
        disp('==========================');

        % save predictions/results/parameters to disk
        filename = ['predscores_' predictionMethod '_' cv_setting '_' int2str(ds) '_' int2str(k) '.mat'];
        save(filename, 'Y_final')
    end

end