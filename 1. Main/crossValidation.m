function aupr=crossValidation(Y)
%crossValidation runs cross validation experiments
%
% INPUT:
%  Y:                     interaction matrix
%
% OUTPUT:
%  scores:                final evaluation score (AUC)

    % Parameters
    global gridSearchMode

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% cross validation setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global cv_setting m n

    if ~gridSearchMode
        % parameters
        fprintf('m = %i\n',m);
        fprintf('n = %i\n',n);
    end

    % seeds
    seeds = [];

    % fixed seeds (for reproducibility)
    % S1 setting: all drug-target interactions (pair prediction)
    if strcmp(cv_setting, 'S1')
        seeds = [7993  252  5054  8720  7519];  % S1 (PAIR)

    % S2 setting: all drugs (drug prediction)
    elseif strcmp(cv_setting, 'S2')
        seeds = [6381  1535  9727  9129  6802]; % S2 (DRUG)

    % S3 setting: all targets (target prediction)
    elseif strcmp(cv_setting, 'S3')
        seeds = [8768  7952  402  9781  9541];  % S3 (TARGET)

    else
        disp('Warning: CV setting not recognized!')
        disp(' ')
    end

    % if no seeds specified...
    if isempty(seeds)
        rng('shuffle');
        seeds = randi(10000,1,m);
    end

    if ~gridSearchMode
        % print seeds to be used...
        fprintf('seeds = [  ');
        for s=1:length(seeds)
            fprintf('%i  ',seeds(s));
        end
        fprintf('  ]\n\n');
        disp('==========================');
        diary off; diary on;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % cross validation (m repetitions of n-fold experiments)
    AUCs  = zeros(1,m);
    AUPRs = zeros(1,m);
    for i=1:m
        seed = seeds(i);
        [AUCs(i), AUPRs(i)] = nfold(Y,seed,i);
    end

    if ~gridSearchMode
        % display evaluation results
        fprintf('\n FINAL AVERAGED RESULTS\n\n');
        fprintf('     AUC (std): %g\t(%g)\n',   mean(AUCs),  std(AUCs));
        fprintf('    AUPR (std): %g\t(%g)\n',   mean(AUPRs), std(AUPRs));
    end

    aupr = mean(AUPRs);
end