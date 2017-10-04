function [patterns,labels,train_indices]=generateTrainingSet(Y,exclude_indices,drugFeatVectors,targetFeatVectors)
%generateTrainingSet generates a training set to be used for training
%classifiers
%
% INPUT:
%  Y:                   interaction matrix
%  exclude_indices:     indices of instances to be excluded from training set
%  drugFeatVectors:     drug feature vectors
%  targetFeatVectors:   target feature vectors
%
% OUTPUT:
%  patterns:            feature vectors of the training set instances
%  labels:              corresponding labels of the training set instances
%  train_indices:       indices of training set instances

    % +ve set
    pos_indices = find(Y == 1);                             % +ve instances
    pos_indices = setdiff(pos_indices, exclude_indices);    % remove instances to be excluded (e.g. test set instances)
    pos_indices = pos_indices(:);                           % force into column vector

    % -ve set
    global npRatio                                                                                  % -ve to +ve ratio (specified in 'begin.m' or 'begin_multiple.m')
    neg_indices = find(Y ~= 1);                                                                     % -ve instances
    neg_indices = setdiff(neg_indices, exclude_indices);                                            % remove instances to be excluded (e.g. test set instances)
    neg_indices = neg_indices(randperm(length(neg_indices), ceil(npRatio*length(pos_indices))));    % random undersampling for -ve set
    neg_indices = neg_indices(:);                                                                   % force into column vector

    % training instances: generate feature vectors, labels and weights
    posFeatureVectors = generateFeatures(pos_indices, drugFeatVectors, targetFeatVectors);  % +ve feature vectors
    negFeatureVectors = generateFeatures(neg_indices, drugFeatVectors, targetFeatVectors);  % -ve feature vectors
    patterns = [posFeatureVectors; negFeatureVectors];                                      % training feature vectors
    labels   = [ones(length(pos_indices),1); zeros(length(neg_indices),1)];                 % training labels
    train_indices = [pos_indices; neg_indices];                                             % indices of training instances
end