%-------------------------------------------------------------------
% Before running this file:
% Run an experiment using start.m. This will lead to .mat file being saved 
% containing the results. After doing so, load the .mat file and THEN run
% THIS file to show the previously computed results
% 
% e.g. run start.m with 'tree' as the prediction method for a quick trial.
% 
% NOTE:
% For more details, on what is stored in the .mat file, check
% crossValidation.m
%-------------------------------------------------------------------

% load data
Y = importdata('data/interactionMatrix.txt');


% seeds used to generate results in paper
cv_setting = 'S1';
switch cv_setting
    case 'S1',  seed = 13579;  % S1 (PAIR)
    case 'S2',  seed = 8148;   % S2 (DRUG)
    case 'S3',  seed = 9706;   % S3 (TARGET)
end


N = 5;                                              % number of folds
CV_indices = getCVindices(seed,cv_setting,N,Y);     % get CV indices


% to hold evaluation scores for each of the N folds
aucs = zeros(N,1);

% load results
for n=1:N
    fprintf('\nFOLD %g:\t',n);

    % test set
    test_indices = find(CV_indices==n);

    % evaluate predictions
    aucs(n) = calculate_auc(predScoresMatrix(test_indices), Y(test_indices));
    fprintf('%.3g', aucs(n));
end

fprintf('\n\nAUC:\t%.3g\n\n', mean(aucs));
diary off;

%-------------------------------------------------------------------