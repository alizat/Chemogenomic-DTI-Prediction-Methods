clear all

%-------------------------------------------------------------------

warning off
diary off; diary on;
fprintf('\nSTART TIME:    %s\n\n', datestr(now));

%-------------------------------------------------------------------

global drugFeatureVectors targetFeatureVectors

% load data
path = 'data/';
Y = importdata([path 'interactionMatrix.txt']);
drugFeatureVectors = importdata([path 'drugFeatureVectors.txt']);
targetFeatureVectors = importdata([path 'targetFeatureVectors.txt']);

%-------------------------------------------------------------------

global batchSize cv_setting gridSearchMode

batchSize = 10000;      % batch size
gridSearchMode = 0;     % grid search mode?

%-------------------------------------------------------------------

global predictionMethod npRatio numLearners r

predictionMethods = {'dt','rf','ensemdt'};    % ,'svm'  (NOTE: SVM takes quite a long time to execute!)
for pm=1:length(predictionMethods)
    disp('==================================================')
    predictionMethod = predictionMethods{pm};
    if ismember(predictionMethod, {'dt', 'svm', 'rf'})
        npRatio = 1;            % -ve to +ve ratio
        
    elseif strcmp(predictionMethod, 'ensemdt')
        npRatio = 5;            % -ve to +ve ratio
        numLearners = 50;       % number of base learners
        r = 0.2;                % percentage of features to be used per base learner (feature subspacing)
    else
        error(['Unrecognized method: '  predictionMethod])
    end

    % loop over all cross validation settings
    for cvs=1:3
        disp('-------------------')

        % CV setting
        cv_setting = ['S' int2str(cvs)];
        switch cv_setting
            case 'S1', disp('CV Setting Used: S1 - PAIR');
            case 'S2', disp('CV Setting Used: S2 - DRUG');
            case 'S3', disp('CV Setting Used: S3 - TARGET');
        end
        disp(' ')

        %-----------------------------

        % print parameters
        %disp(['     batchSize = ' num2str(batchSize)])
        %disp(['gridSearchMode = ' num2str(gridSearchMode)])
        %disp(' ')
        disp(['Prediction method = ' predictionMethod])
        disp(['          npRatio = ' num2str(npRatio)])
        if strcmp(predictionMethod, 'ensemble')
            disp(['      numLearners = ' num2str(numLearners)])
            disp(['                r = ' num2str(r)])
        end

        %-----------------------------

        % run chosen selection method and output CV results
        tic
        scores = crossValidation(Y);
        disp(' ')
        toc
        fprintf('\n\nAUC:\t%.3g\n\n\n', scores.auc)

        disp('-------------------')
        disp(' ')
        disp(' ')
        diary off; diary on;
    end
    
    disp('==================================================')
    diary off; diary on;
end
diary off;

%-------------------------------------------------------------------