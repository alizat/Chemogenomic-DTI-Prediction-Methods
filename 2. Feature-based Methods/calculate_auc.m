function auc=calculate_auc(predictionScores,labels)
    labels = (labels > 0) - (labels <= 0);
    [~,idx] = sort(predictionScores, 'descend');
    labels = labels(idx);
    tp = cumsum(labels == 1);
    fp = sum(labels == -1);
    auc = sum(tp(labels == -1));
    if tp == 0 | fp == 0
        %fprintf('warning: Too few postive true labels or negative true labels')
        auc = 0;
    else
        auc = auc / tp(end) / fp;
    end
end


% function auc = calculate_auc(targets, predicts)
% 	% Calculate area under the ROC curve
% 	if nargin > 1
% 		[~,i] = sort(predicts,'descend');
% 		targets = targets(i);
% 	end
% 	
% 	%for i=1:n
% 	%	if targets(i)
% 	%		goods = goods + 1
% 	%	else
% 	%		auc = auc + goods;
% 	%	end
% 	%end
% 	cumsums = cumsum(targets);
% 	auc = sum(cumsums(~targets));
% 	pos = sum(targets);
% 	neg = sum(~targets);
% 	if pos == 0, warning('Calculate auc: no positive targets'); end
% 	if neg == 0, warning('Calculate auc: no negative targets'); end
% 	auc = auc / (pos * neg + eps);
% end