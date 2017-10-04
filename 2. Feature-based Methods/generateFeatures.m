function featureVectors=generateFeatures(instances,drugFeatVectors,targetFeatVectors)
%generateFeatures generates the instance feature vectors for the supplied
%indices
%
% INPUT:
%  instaces:            interaction matrix
%  drugFeatVectors      feature vectors of drugs
%  targetFeatVectors    feature vectors of targets
% 
% OUTPUT:
%  featureVectors:      instance feature vectors

    % prepare drug-target pairs
    numDrugs = size(drugFeatVectors, 1);
    numTargets = size(targetFeatVectors, 1);
    [d,t] = ind2sub([numDrugs numTargets], instances);

    % form concatenated feature vectors
    featureVectors = drugFeatVectors(d, :);                      % drug columns
    featureVectors = [featureVectors targetFeatVectors(t, :)];   % target columns

%     % equivalent of the above two lines, but more readable (and slower)
%     numDrugFeatures = size(drugFeatVectors, 2);
%     numTargetFeatures = size(targetFeatVectors, 2);
%     featureVectors = zeros(length(instances), numDrugFeatures + numTargetFeatures);
%     for i=1:size(instances, 1)
%         featureVectors(i, :) = [drugFeatVectors(d(i), :)  targetFeatVectors(t(i), :)];
%     end

end