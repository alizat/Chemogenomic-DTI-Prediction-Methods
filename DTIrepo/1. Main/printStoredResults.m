%-------------------------------------------------------------------
% Before running this file:
% Run an experiment using start.m. This will lead to .mat file being saved 
% containing the results. After doing so, load the .mat file and THEN run
% THIS file to show the previously computed results
% 
% e.g. run start.m with 'np' as the prediction method for a quick trial.
% 
% NOTE:
% For more details on what is stored in the .mat file, check
% crossValidation.m
%-------------------------------------------------------------------

% load data
Y = importdata('data/interactionMatrix.txt');

% load saved matlab files
load('CV_indices.mat')

% number of folds
N = length(unique(CV_indices));

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