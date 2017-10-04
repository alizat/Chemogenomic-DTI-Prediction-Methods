function Yhat=alg_template(Y,predictionMethod,test_indices)
%alg_template predicts DTIs based on the prediction method selected in
%start.m or sensitivity_analysis.m
%
% INPUT:
%  Y:                     interaction matrix
%  predictionMethod:      method to use for prediction
%  test_indices:          indices of the test set instances
%
% OUTPUT:
%  Yhat:                  prediction scores matrix

    % Parameters
    global drugFeatureVectors targetFeatureVectors

    switch predictionMethod
        % proposed ensemble framework --------------------
        case 'ensemdt'
            Yhat = alg_ensemble(Y, test_indices);

        % rest of feature-based methods ------------------
        otherwise
            % generate training set
            [patterns, labels, ~] = generateTrainingSet(Y, test_indices, drugFeatureVectors, targetFeatureVectors);

            % train predictive model
            switch predictionMethod
                case 'dt',      predModel.model = compact(fitctree(patterns, labels, 'MinLeafSize', 10));
                case 'rf',      predModel.model = compact(TreeBagger(500, patterns, labels, 'Prior', 'uniform'));
                case 'svm',     predModel.model = compact(fitcsvm(patterns, labels, 'KernelFunction', 'rbf'));
            end
            clear patterns labels

            % predict
            Yhat = predictor(predModel, test_indices, predictionMethod, drugFeatureVectors, targetFeatureVectors);
    end

end


% *********************************************************************** %
% *********************************************************************** %
% *********************************************************************** %
% *********************************************************************** %
% *********************************************************************** %
% *********************************************************************** %
% *********************************************************************** %
% *********************************************************************** %
% *********************************************************************** %
% *********************************************************************** %


function Yhat=alg_ensemble(Y,test_indices)
%alg_ensemble predicts DTIs based on the algorithm described in the
%following paper (but without the dimensionality reduction):
% Ali Ezzat, Peilin Zhao, Min Wu, Xiao-Li Li and Chee-Keong Kwoh
% (2017) Drug-Target Interaction Prediction using Ensemble Learning and
%           Dimensionality Reduction 
%
% INPUT:
%  Y:                     interaction matrix
%  test_indices:          indices of the test set instances
%
% OUTPUT:
%  Yhat:                  prediction matrix (including test set scores)

    % Parameters
    global numLearners drugFeatureVectors targetFeatureVectors r

    % instances to be excluded from training set
    exclude_indices = test_indices;

    % generate models
    Yhat = 0;
    for c = 1:numLearners
        % feature subspacing
        numDrugFeatures = size(drugFeatureVectors, 2);
        numTargetFeatures = size(targetFeatureVectors, 2);
        drugFeatures = randperm(numDrugFeatures, floor(numDrugFeatures*r));
        targetFeatures = randperm(numTargetFeatures, floor(numTargetFeatures*r));
        drugFeatVectors = drugFeatureVectors(:, drugFeatures);
        targetFeatVectors = targetFeatureVectors(:, targetFeatures);

        % train
        [patterns, labels, ~] = generateTrainingSet(Y, exclude_indices, drugFeatVectors, targetFeatVectors);
        predModel.model = compact(fitctree(patterns, labels, 'Prior', 'uniform'));
        Yhat = Yhat + predictor(predModel, test_indices, 'ensemdt', drugFeatVectors, targetFeatVectors);
    end
    clear patterns labels drugFeatVectors targetFeatVectors
end


% *********************************************************************** %
% *********************************************************************** %
% *********************************************************************** %
% *********************************************************************** %
% *********************************************************************** %
% *********************************************************************** %
% *********************************************************************** %
% *********************************************************************** %
% *********************************************************************** %
% *********************************************************************** %


function Yhat=predictor(predModel,test_indices,algorithm,drugFeatVectors,targetFeatVectors)
%predictor is a helper function that produces prediction scores for
%feature-based DTI prediction methods

    global batchSize

    % to hold prediction scores
    predScores = zeros(length(test_indices), 1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTE:                                                            %
    % the entire set of feature vectors of the test set can't fit into %
    % memory, so we do predictions of test set instances in batches    %
    % (default batch size = 10000)                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % collect predictions
    for counter = 0:batchSize:length(test_indices)
        % current batch of test set instances
        start_index = counter + 1;
        end_index = min(counter + batchSize, length(test_indices));
        if start_index > length(test_indices), break; end

        % prepare testing instances
        testingFeatureVectors = generateFeatures(test_indices(start_index:end_index), drugFeatVectors, targetFeatVectors);

        % get prediction scores
        predScores(start_index:end_index) = predictWithModel(algorithm, predModel, testingFeatureVectors);

        clear testingFeatureVectors
    end

    Yhat = zeros(size(drugFeatVectors, 1), size(targetFeatVectors, 1));
    Yhat(test_indices) = predScores;
end

% ----------------------------------------------------------------------

function scores=predictWithModel(algorithm,predModel,X)

    switch algorithm
        % DECISION TREE, RANDOM FOREST, SVM ----------
        case {'dt', 'rf', 'svm'}
            % prediction scores (probabilities)
            [~, scores] = predict(predModel.model, X);
            scores = scores(:,2);

        % OUR METHOD ---------------------------------
        otherwise
            % predicted labels
            scores = predict(predModel.model, X);
            if iscell(scores), scores = str2double(scores); end
    end
end