function scores=crossValidation(Y)
%crossValidation runs cross validation experiments
%
% INPUT:
%  Y:                     interaction matrix
%
% OUTPUT:
%  scores:                final evaluation score (AUC)

    % Parameters
    global predictionMethod Yn

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% cross validation setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global cv_setting gridSearchMode
    N = 5;  % number of folds in the N-fold CV

    if ~gridSearchMode
        % parameters
        fprintf('N = %i\n',N);
    end

    % seeds
    seed = [];

    % seeds used to generate results
    % S1 setting: all drug-target interactions (pair prediction)
    if strcmp(cv_setting, 'S1')
        seed = 13579;  % S1 (PAIR)

    % S2 setting: all drugs (drug prediction)
    elseif strcmp(cv_setting, 'S2')
        seed = 8148;   % S2 (DRUG)

    % S3 setting: all targets (target prediction)
    elseif strcmp(cv_setting, 'S3')
        seed = 9706;   % S3 (TARGET)
    end

    % if no seeds specified...
    if isempty(seed)
        rng('shuffle');
        seed = randi(10000);
    end

    if ~gridSearchMode
        % print seeds to be used...
        fprintf('seeds = [  ');
        for s=1:length(seed)
            fprintf('%i  ',seed(s));
        end
        fprintf('  ]\n\n');
        disp('==========================');
        diary off; diary on;
    end

    % get CV indices
    CV_indices = getCVindices(seed,cv_setting,N,Y);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % to pool predictions from all N folds
    predScoresMatrix = zeros(size(Y));  % predicted scores (continuous)

    % to hold predictions for test indices in each of the N folds
    predScores = cell(N,1);

    % to hold evaluation scores for each of the N folds
    aucs = zeros(N,1);

    diary off; diary on;
    for n=1:N
        if ~gridSearchMode
            fprintf('\nFOLD %g:\t',n);
        end

        % test set
        test_indices = find(CV_indices==n);
        test_indices = test_indices(:);

        % leave out interactions
        Yn = Y;
        Yn(test_indices) = 0;

        % to ensure reproducibility of results 
        % (random generation of -ve's in generateNegativeSet.m)
        rng(n * 1111)

        % predict
        Yhat = alg_template(Yn, predictionMethod, test_indices);
        predScores{n} = Yhat(test_indices);
        predScoresMatrix(test_indices) = predScores{n};

        % evaluate predictions
        aucs(n) = calculate_auc(predScores{n}, Y(test_indices));

        % print results
        if ~gridSearchMode
            fprintf('%.3g\t\t\t\tTIME:    %s', aucs(n), datestr(now));
        end
        
        diary off; diary on;

        if gridSearchMode
            aucs = aucs(n);
            break
        end
    end

    if ~gridSearchMode
        % save prediction results
        currentTime = strrep(datestr(now), ' ', '_');
        currentTime = strrep(currentTime, ':', '_');
        filename = ['predscores_' predictionMethod '_' cv_setting '_' currentTime '.mat'];

        % save predictions/results/parameters to disk
        save(filename, 'predScoresMatrix')
    end

    % final averaged scores
    scores.auc = mean(aucs);
end