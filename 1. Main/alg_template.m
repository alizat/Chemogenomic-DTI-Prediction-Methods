function Yhat=alg_template(Y,predictionMethod,test_ind,left_out)
%alg_template predicts DTIs based on the prediction method selected in
%start.m or sensitivity_analysis.m
%
% INPUT:
%  Y:                     interaction matrix
%  predictionMethod:      method to use for prediction
%  test_indices:          indices of the test set instances
%  left_out:              in case of S1: left_out = test_indices
%                         in case of S2: left_out = left out drugs
%                         in case of S3: left_out = left out targets
%
% OUTPUT:
%  Yhat:                  prediction scores matrix

    % Parameters
    global Sd St
    predFn = str2func(['alg_'  predictionMethod]);
    Yhat = predFn(Y,Sd,St,test_ind,left_out);

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


function Yhat=alg_np(Y,Sd,St,~,~)
%alg_np predicts DTIs based on the Nearest Profile algorithm described in the following paper: 
% Yoshihiro Yamanishi, Michihiro Araki, Alex Gutteridge, Wataru Honda and Minoru Kanehisa,
% (2008) Prediction of drug–target interaction networks from the integration of chemical and genomic spaces

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Nearest Profile (NP) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Drug
    Sd(logical(eye(length(Sd)))) = 0;   % remove self-similarities
    [maxx, indx] = max(Sd);             % get nearest neighbor for each drug
    for i=1:length(Sd)
        Sd(i, :) = 0;                   % reset all similarities to 0...
        Sd(i, indx(i)) = maxx(i);       % except that of the nearest neighbor
    end
    yd = Sd * Y;

    % Target
    St(logical(eye(length(St)))) = 0;   % remove self-similarities
    [maxx, indx] = max(St);             % get nearest neighbor for each target
    for j=1:length(St)
        St(j, :) = 0;                   % reset all similarities to 0...
        St(j, indx(j)) = maxx(j);       % except that of the nearest neighbor
    end
    yt = (St * Y')';

    % Final Result
    Yhat = (yd + yt) / 2;
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


function Yhat=alg_wp(Y,Sd,St,~,~)
%alg_wp predicts DTIs based on the Weighted Profile algorithm described in the following paper: 
% Yoshihiro Yamanishi, Michihiro Araki, Alex Gutteridge, Wataru Honda and Minoru Kanehisa,
% (2008) Prediction of drug–target interaction networks from the integration of chemical and genomic spaces

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Weighted Profile (WP) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    yd = bsxfun(@rdivide, Sd * Y, sum(Sd,2));   yd(Y==1) = 1;   % Drug
    yt = bsxfun(@rdivide, Y * St, sum(St));     yt(Y==1) = 1;   % Target
    Yhat = (yd + yt) / 2;
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


function Yhat=alg_rls_wnn(Y,ka,kb,~,~)
%alg_rls_wnn predicts DTIs based on the algorithm described in the following paper: 
% Twan van Laarhoven, Elena Marchiori,
% (2013) Predicting drug–target interactions for new drug compounds using a
%           weighted nearest neighbor profile 
% 
% Code below is adapted from the code available at this website:
% http://cs.ru.nl/~tvanlaarhoven/drugtarget2013/

    %%%%%%%%%%%
    %%% WNN %%%
    %%%%%%%%%%%

    global eta
    %eta = 0.7;     %default
    Y = preprocess_WNN(Y,ka,kb,eta);


    %%%%%%%%%%%
    %%% GIP %%%
    %%%%%%%%%%%

    global alpha
    %alpha = 0.5;   %default
    ka = alpha*ka + (1-alpha)*getGipKernel(Y);
    kb = alpha*kb + (1-alpha)*getGipKernel(Y');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Regularized Least Squares (RLS-kron) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global sigma
	%sigma = 1;     %default
	[va,la] = eig(ka);
	[vb,lb] = eig(kb);
	l = kron(diag(lb)',diag(la));
	l = l ./ (l + sigma);
	m1 = va' * Y * vb;
	m2 = m1 .* l;
	Yhat = va * m2 * vb';
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


function Yhat=alg_nbi(Y,Sd,St,~,~)
%alg_nbi predicts DTIs based on the algorithm described in the following paper: 
% Feixiong Cheng, Chuang Liu, Jing Jiang, Weiqiang Lu, Weihua Li, Guixia Liu, Weixing Zhou, Jin Huang, Yun Tang
% (2012) Prediction of Drug-Target Interactions and Drug Repositioning via Network-Based Inference

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Network-based Inference (NBI) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % normalize Sd and St:
    Sd = Sd ./ (sum(Sd,2) * sum(Sd));
    St = St ./ (sum(St,2) * sum(St));
    % based on Equation (3) from:
    % Wenhui Wang, Sen Yang, Jing Li
    % (2013) Drug target predictions based on heterogeneous graph inference

    % NBI
    global alpha
    %alpha = 0.5;   %default
    Yhat = Y;
    Yhat = (alpha * Sd * Yhat * St) + ((1 - alpha) * Y);
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


function y3=alg_kbmf2k(Y,Sd,St,~,~)
%alg_kbmf2k predicts DTIs based on the algorithm described in the following paper:
% Mehmet Gönen
% (2012) Predicting drug–target interactions from chemical and genomic kernels using Bayesian matrix factorization

    %--------------------------------------------------------------------

    % parameters
    global rs
    params.R = rs;

    addpath('kbmf2k');

    %--------------------------------------------------------------------
	
	%  where yij = +1 if drug compound di interacts with target protein tj and yij = -1 otherwise
	Y = 2 * (Y > 0.5) - 1;

    % training
    state = alg_kbmf_regression_train(Sd, St, Y, params);

    % predict
    y3 = alg_kbmf_regression_test(Sd, St, state);

end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Mehmet Gonen (mehmet.gonen@aalto.fi)
% Helsinki Institute for Information Technology HIIT
% Department of Information and Computer Science
% Aalto University School of Science

function prediction = alg_kbmf_regression_test(Kx, Kz, state)
%     addpath('kbmf2k');
%     addpath('kbmf1mkl1k');
%     addpath('kbmf1k1mkl');
%     addpath('kbmf2mkl');
    prediction = state.parameters.test_function(Kx, Kz, state);
    prediction = prediction.Y.mean;
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Mehmet Gonen (mehmet.gonen@aalto.fi)
% Helsinki Institute for Information Technology HIIT
% Department of Information and Computer Science
% Aalto University School of Science

function state = alg_kbmf_regression_train(Kx, Kz, Y, otherParameters, varargin)
%     addpath('kbmf2k');
%     addpath('kbmf1mkl1k');
%     addpath('kbmf1k1mkl');
%     addpath('kbmf2mkl');


    % external parameters
    if isfield(otherParameters,'R'), parameters.R = otherParameters.R; end
    if isfield(otherParameters,'margin'), parameters.margin = otherParameters.margin; end


    Px = size(Kx, 3);
    Pz = size(Kz, 3);
    %is_supervised = all(~isnan(Y(:)));
    is_supervised = 0;

    parameters.alpha_lambda = 1;
    parameters.beta_lambda = 1;
    if Px > 1 || Pz > 1
        parameters.alpha_eta = 1;
        parameters.beta_eta = 1;
    end
    parameters.iteration = 50;
    parameters.progress = 1;
    parameters.seed = 1606;
    parameters.sigmag = 0.1;
    if Px > 1 || Pz > 1
        parameters.sigmah = 0.1;
    end
    parameters.sigmay = 1.0;

    if is_supervised == 1
        if Px == 1 && Pz == 1
            train_function = @kbmf2k_supervised_regression_variational_train_with_bound;
            test_function = @kbmf2k_supervised_regression_variational_test;
        elseif Px > 1 && Pz == 1
            train_function = @kbmf1mkl1k_supervised_regression_variational_train;
            test_function = @kbmf1mkl1k_supervised_regression_variational_test;
        elseif Px == 1 && Pz > 1
            train_function = @kbmf1k1mkl_supervised_regression_variational_train;
            test_function = @kbmf1k1mkl_supervised_regression_variational_train;
        elseif Px > 1 && Pz > 1
            train_function = @kbmf2mkl_supervised_regression_variational_train;
            test_function = @kbmf2mkl_supervised_regression_variational_test;
        end
    else
        if Px == 1 && Pz == 1
            train_function = @kbmf2k_semisupervised_regression_variational_train;
            test_function = @kbmf2k_semisupervised_regression_variational_test;
        elseif Px > 1 && Pz == 1
            train_function = @kbmf1mkl1k_semisupervised_regression_variational_train;
            test_function = @kbmf1mkl1k_semisupervised_regression_variational_test;
        elseif Px == 1 && Pz > 1
            train_function = @kbmf1k1mkl_semisupervised_regression_variational_train;
            test_function = @kbmf1k1mkl_semisupervised_regression_variational_test;
        elseif Px > 1 && Pz > 1
            train_function = @kbmf2mkl_semisupervised_regression_variational_train;
            test_function = @kbmf2mkl_semisupervised_regression_variational_test;
        end
    end

    for i = 1:2:nargin - 4
        parameters.(varargin{i}) = varargin{i + 1};
    end
    
    parameters.train_function = train_function;
    parameters.test_function = test_function;

    state = train_function(Kx, Kz, Y, parameters);
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


function y3=alg_grmf(Y,Sd,St,test_ind,~)
%alg_grmf predicts DTIs based on the algorithm described in the following paper: 
% Ali Ezzat, Peilin Zhao, Min Wu, Xiao-Li Li and Chee-Keong Kwoh
% (2016) Drug-target interaction prediction with graph-regularized matrix factorization
%
% INPUT:
%  Y:           interaction matrix
%  Sd:          pairwise row similarities matrix
%  St:          pairwise column similarities matrix
%  test_ind:    indices of test set instances
%
% OUTPUT:
%  y3:          prediction matrix
%

    % parameters
    global num_iter p k lambda_l lambda_d lambda_t

    % preprocessing Sd & St
    % (Sparsification of matrices via p-nearest-neighbor graphs)
    Sd = preprocess_PNN(Sd,p);
    St = preprocess_PNN(St,p);

    % Laplacian Matrices
    Dd = diag(sum(Sd));
    Ld = Dd - Sd;
    Ld = (Dd^(-0.5))*Ld*(Dd^(-0.5));
    Dt = diag(sum(St));
    Lt = Dt - St;
    Lt = (Dt^(-0.5))*Lt*(Dt^(-0.5));

    % (W)GRMF
    [A,B] = initializer(Y,k);	% initialize A & B
    W = ones(size(Y));          % weight matrix W
    W(test_ind) = 0;            % set W=0 for test instances
    [A,B] = alg_grmf_predict(Y,A,B,Ld,Lt,lambda_l,lambda_d,lambda_t,num_iter,W);    % update A & B

    % compute prediction matrix
    y3 = A*B';

end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [A,B]=alg_grmf_predict(Y,A,B,Ld,Lt,lambda_l,lambda_d,lambda_t,num_iter,W)
%alg_grmf_predict performs alternating least squares for GRMF
%
% INPUT:
%  Y:           interaction matrix
%  A:           drug latent feature matrix
%  B:           target latent feature matrix
%  Ld:          drug graph Laplacian
%  Lt:          target graph Laplacian
%  lambda_ldt:  regularization parameters
%  num_iter:    number of iterations for alternating least squares
%  W:           weight matrix
%
% OUTPUT:
%  A:           updated drug latent feature matrix
%  B:           updated target latent feature matrix
%
    
    K = size(A,2);
    lambda_d_Ld = lambda_d*Ld;          % to avoid 
    lambda_t_Lt = lambda_t*Lt;          % repeated matrix 
    lambda_l_eye_K = lambda_l*eye(K);   % multiplications

    % if no weight matrix is supplied or W is an all-ones matrix...
    if nargin < 10 || isequal(W,ones(size(W)))

        %%%%%%%%%%%%
        %%% GRMF %%%
        %%%%%%%%%%%%

        for z=1:num_iter
            A = (Y*B  - lambda_d_Ld*A) / (B'*B + lambda_l_eye_K);
            B = (Y'*A - lambda_t_Lt*B) / (A'*A + lambda_l_eye_K);
        end
        
        
    else

        %%%%%%%%%%%%%
        %%% WGRMF %%%
        %%%%%%%%%%%%%

        H = W .* Y;
        for z=1:num_iter
%             % for readability...
%             A_old = A;
%             for i=1:size(A,1)
%                 A(i,:) = (H(i,:)*B - lambda_d*Ld(i,:)*A_old) / (B'*diag(W(i,:))*B + lambda*eye(k));
%             end
%             B_old = B;
%             for j=1:size(B,1)
%                 B(j,:) = (H(:,j)'*A - lambda_t*Lt(j,:)*B_old) / (A'*diag(W(:,j))*A + lambda*eye(k));
%             end

            % equivalent, less readable, faster
            A_old = A;
            HB_minus_alpha_Ld_A_old = H*B - lambda_d_Ld*A_old;
            for a=1:size(A,1)
                A(a,:) = HB_minus_alpha_Ld_A_old(a,:) / (B'*diag(W(a,:))*B + lambda_l_eye_K);
            end
            
            B_old = B;
            HtA_minus_beta_Lt_B_old = H'*A - lambda_t_Lt*B_old;
            for b=1:size(B,1)
                B(b,:) = HtA_minus_beta_Lt_B_old(b,:) / (A'*diag(W(:,b))*A + lambda_l_eye_K);
            end
        end
    end
    
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [S,p]=preprocess_PNN(S,p)
%preprocess_PNN sparsifies S by keeping, for each drug/target, the "p"
% nearest neighbors (NNs) and discarding the rest. 

    NN_mat = zeros(size(S));
    for j=1:length(NN_mat)
        [~,indx] = sort(S(j,:),'descend');
        indx = indx(1:p+1);     % keep drug/target j and its "p" NNs
        NN_mat(j,indx) = 1;
    end
    NN_mat = (NN_mat+NN_mat')/2;
    S = NN_mat .* S;

end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [A,B]=initializer(Y,k)
%initializer initializes the A and B latent feature matrices for either
% of the CMF or GRMF algorithms.
%
% INPUT:
%  Y:   interaction matrix
%  k:   number of latent features
%
% OUTPUT:
%  A:   latent feature matrix for drugs
%  B:   latent feature matrix for targets
%

    [u,s,v] = svds(Y,k);
    A = u*(s^0.5);
    B = v*(s^0.5);

%     % Alternative: Use non-negative matrix factorization
%     k = min(k, min(size(Y)));
%     [A,B] = nnmf(Y,k);
%     B = B';

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

function y3=alg_cmf(Y,Sd,St,test_ind,~)
%alg_cmf predicts DTIs based on the algorithm described in the following paper:
% Xiaodong Zheng, Hao Ding, Hiroshi Mamitsuka and Shanfeng Zhu
% (2013) Collaborative Matrix Factorization with Multiple Similarities for Predicting Drug-Target Interactions
%
% INPUT:
%  Y:           interaction matrix
%  Sd:          pairwise row similarities matrix
%  St:          pairwise column similarities matrix
%  test_ind:    indices of test set instances
%
% OUTPUT:
%  y3:          prediction matrix
%

    % parameters
    global num_iter k lambda_l lambda_d lambda_t

    [A,B] = initializer(Y,k);	% initialize A & B
    W = ones(size(Y));          % weight matrix W
    W(test_ind) = 0;            % set W=0 for test instances
    [A,B] = alg_cmf_predict(Y,A,B,Sd,St,lambda_l,lambda_d,lambda_t,num_iter,W);     % update A & B

    % compute prediction matrix
    y3 = A*B';

end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [A,B]=alg_cmf_predict(Y,A,B,Sd,St,lambda_l,lambda_d,lambda_t,num_iter,W)
%alg_cmf_predict performs alternating least squares for CMF
%
% INPUT:
%  Y:           interaction matrix
%  A:           drug latent feature matrix
%  B:           target latent feature matrix
%  Sd:          pairwise row similarities matrix
%  St:          pairwise column similarities matrix
%  lambda_ldt:  regularization parameters
%  num_iter:    number of iterations for alternating least squares
%  W:           weight matrix
%
% OUTPUT:
%  A:           updated drug latent feature matrix
%  B:           updated target latent feature matrix
%
    
    K = size(A,2);
    lambda_d_Sd = lambda_d*Sd;          % to avoid 
    lambda_t_St = lambda_t*St;          % repeated matrix 
    lambda_l_eye_K = lambda_l*eye(K);   % multiplications

    % if no weight matrix is supplied or W is an all-ones matrix...
    if nargin < 10 || isequal(W,ones(size(W)))
        AtA = A'*A;
        BtB = B'*B;
        for z=1:num_iter
            A = (Y*B + lambda_d_Sd*A)  / (BtB + lambda_l_eye_K + lambda_d*(AtA));
            AtA = A'*A;
            B = (Y'*A + lambda_t_St*B) / (AtA + lambda_l_eye_K + lambda_t*(BtB));
            BtB = B'*B;
        end
        
    else
        H = W .* Y;
        for z=1:num_iter
%             % for readability...
%             A_old = A;
%             lambda_d_A_oldt_A_old = lambda_d*(A_old'*A_old);
%             for a=1:size(A,1)
%                 A(a,:) = (H(a,:)*B + lambda_d_Sd(a,:)*A_old) / (B'*B + lambda_l_eye_k + lambda_d_A_oldt_A_old);
%             end
%             B_old = B;
%             lambda_t_B_oldt_B_old = lambda_t*(B_old'*B_old);
%             for b=1:size(B,1)
%                 B(b,:) = (H(:,b)'*A + lambda_t_St(b,:)*B_old) / (A'*A + lambda_l_eye_k + lambda_t_B_oldt_B_old);
%             end

            % equivalent, less readable, faster
            A_old = A;
            HB_plus_lambda_d_Sd_A_old = H*B + lambda_d_Sd*A_old;
            lambda_l_eye_k_plus_lambda_d_A_oldt_A_old = lambda_l_eye_K + lambda_d*(A_old'*A_old);
            for a=1:size(A,1)
                A(a,:) = HB_plus_lambda_d_Sd_A_old(a,:) / (B'*diag(W(a,:))*B + lambda_l_eye_k_plus_lambda_d_A_oldt_A_old);
            end
            B_old = B;
            HtA_plus_lambda_t_St_B_old = H'*A + lambda_t_St*B_old;
            lambda_l_eye_k_plus_lambda_t_B_oldt_B_old = lambda_l_eye_K + lambda_t*(B_old'*B_old);
            for b=1:size(B,1)
                B(b,:) = HtA_plus_lambda_t_St_B_old(b,:) / (A'*diag(W(:,b))*A + lambda_l_eye_k_plus_lambda_t_B_oldt_B_old);
            end

        end
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


function Yhat=alg_sitar(Y,Sd,St,test_ind,~)
%alg_nbi predicts DTIs based on the algorithm described in the following paper: 
% Feixiong Cheng, Chuang Liu, Jing Jiang, Weiqiang Lu, Weihua Li, Guixia Liu, Weixing Zhou, Jin Huang, Yun Tang
% (2012) Prediction of Drug-Target Interactions and Drug Repositioning via Network-Based Inference
%
% INPUT:
%  Y:           interaction matrix
%  Sd:          pairwise row similarities matrix
%  St:          pairwise column similarities matrix
%  test_ind:    indices of test set instances
%
% OUTPUT:
%  Yhat:        prediction matrix
%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Similarity-based Inference of drug-TARgets (SITAR) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % +ve set
    pos_indices = find(Y == 1);     % +ve instances
    pos_indices = pos_indices(:);   % force into column vector

    % generate feature vectors
    r = 0.5;
    [ds, ts] = ind2sub(size(Y), 1:numel(Y));
    [pds, pts] = ind2sub(size(Y), pos_indices);
    feat_vectors = (Sd(ds, pds).^r) .* (St(ts, pts).^(1-r));
%     feat_vectors = zeros(numel(Y), length(pos_indices));
%     for i=1:numel(Y)
%         d = ds(i);
%         t = ts(i);
%         feat_vectors(i,:) = (Sd(d, pds).^r) .* (St(t, pts).^(1-r));
%     end

    % train
    train_ind = 1:numel(Y);
    train_ind(test_ind) = [];
    train_ind = train_ind(:);
    model = compact(fitcsvm(feat_vectors(train_ind,:), Y(train_ind), 'KernelFunction', 'rbf'));
    %model = compact(fitcsvm(feat_vectors(train_ind,:), Y(train_ind), 'KernelFunction', 'rbf', 'BoxConstraint', 10));
    %model = compact(TreeBagger(500, feat_vectors(train_ind,:), Y(train_ind)));

    % predict
    [~, scores] = predict(model, feat_vectors(test_ind,:));
    Yhat = Y;
    Yhat(test_ind) = scores(:,2);

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


function Yhat=alg_laprls(Y,Sd,St,~,~)
%alg_laprls predicts DTIs based on the algorithm described in the following paper: 
% Zheng Xia, Ling-Yun Wu, Xiaobo Zhou, Stephen TC Wong,
% (2010) Semi-supervised drug-protein interaction prediction from heterogeneous biological spaces
% 
% Code adapted from supplementary material of Laarhoven 2011
%
% INPUT:
%  Y:           interaction matrix
%  Sd:          pairwise row similarities matrix
%  St:          pairwise column similarities matrix
%
% OUTPUT:
%  Yhat:        prediction matrix
%

    % Parameters as per the above paper
    ga1 = 1;
    gb1 = 1;
    ga2 = 0.01;
    gb2 = 0.01;
    ba  = 0.3;
    bb  = 0.3;

    sa = Sd;
    sb = St;
    ka = Y*Y';
    kb = Y'*Y;
    wa = (ga1*sa + ga2*ka) / (ga1+ga2);
    wb = (gb1*sb + gb2*kb) / (gb1+gb2);

    da = diag(sum(wa));
    db = diag(sum(wb));
    la = da^-0.5 * (da - wa) * da^-0.5;
    lb = db^-0.5 * (db - wb) * db^-0.5;

    fa = wa / (wa + ba*la*wa) * Y;
    fb = wb / (wb + bb*lb*wb) * Y';

    Yhat = (fa+fb') / 2;
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
