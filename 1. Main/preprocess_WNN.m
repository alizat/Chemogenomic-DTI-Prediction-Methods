function Y=preprocess_WNN(Y,Sd,St,eta)
%preprocess_WNN preprocesses the interaction matrix Y using WNN (weighted
%nearest neighbors)
%
% Y = preprocess_WNN(Y,Sd,St,cv_setting,eta)
%
% INPUT:
%  Y:           interaction matrix
%  Sd:          pairwise row similarities matrix
%  St:          pairwise column similarities matrix
%  eta:         decay rate
%
% OUTPUT:
%  Y:           preprocessed interaction matrix
%
%
% Refer to the following paper for more details on WNN:
%   Twan van Laarhoven and Elena Marchiori
%   (2013) Predicting Drug-Target Interactions for New Drug Compounds Using
%           a Weighted Nearest Neighbor Profile 

    global cv_setting
    y = 0;

    % ---------------------------------------

    if ~strcmp(cv_setting,'S3')                 % if not S3 --> infer interaction profiles for empty rows
        yd = Y;
        empty_rows = find(any(Y,2) == 0);       % get indices of empty rows
        w = eta .^ (0:length(Sd)-1);
        w(w < 10^-4) = [];
        k = length(w);
        for r=1:length(empty_rows)
            i = empty_rows(r);                  % i = current empty row
            drug_sim = Sd(i,:);                 % get similarities of drug i to other drugs
            drug_sim(i) = 0;                    % set self-similarity to ZERO
            [~,indx]=sort(drug_sim,'descend');  % sort descendingly
            indx = indx(1:k);
            yd(i,:) = w * Y(indx,:);            % multiply sorted similarities by decreasing decay values
        end
        y = max(yd, y);
    end

    % ---------------------------------------
    
    if ~strcmp(cv_setting,'S2')                     % if not S2 --> infer interaction profiles for empty columns
        yt = Y;
        empty_cols = find(any(Y) == 0);             % get indices of empty columns
        w = eta .^ (0:length(St)-1);
        w(w < 10^-4) = [];
        k = length(w);
        for c=1:length(empty_cols)
            j = empty_cols(c);                      % j = current empty column
            target_sim = St(j,:);                   % get similarities of target j to other targets
            target_sim(j) = 0;                      % set self-similarity to ZERO
            [~,indx]=sort(target_sim,'descend');    % sort descendingly
            indx = indx(1:k);
            yt(:,j) = Y(:,indx) * w';               % multiply sorted similarities by decreasing decay values
        end
        y = max(yt, y);
    end

    % ---------------------------------------
    
    Y = max(Y, y);
    %Y = y;

end