clear all

%-------------------------------------------------------------------

diary off; diary on;
fprintf('\nSTART TIME:    %s\n\n', datestr(now));

%-------------------------------------------------------------------

global predictionMethod gridSearchMode

gridSearchMode = 1;   % grid search mode?
% if turned on:
% - suppresses printing of intermediate fold results on screen
% - prevents saving of prediction scores

warning off

%-------------------------------------------------------------------

global m n Sd St ds cv_setting

% The location of the folder that contains the data
path='data\';

% the different datasets
datasets={'e','ic','gpcr','nr'};

% CV parameters
m = 1;  % number of n-fold experiments (repetitions)
n = 10; % the 'n' in "n-fold experiment"

%-------------------------------------------------------------------

disp(['gridSearchMode = ' num2str(gridSearchMode)])
disp(' ')

diary off; diary on;

%*********************************************************************
global num_iter k lambda_l lambda_d lambda_t p sigma alpha eta rs
%*********************************************************************

%-------------------------------------------------------------------
%-------------------------------------------------------------------
disp(' ')
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp(' ')
%-------------------------------------------------------------------
%-------------------------------------------------------------------

predictionMethod = 'rls_wnn';

for cvs=1:3
    disp('==============================================================');
    disp(['Prediction method = ' predictionMethod])
    cv_setting = ['S' int2str(cvs)];
    switch cv_setting
        case 'S1', disp('CV Setting Used: S1 - PAIR');
        case 'S2', disp('CV Setting Used: S2 - DRUG');
        case 'S3', disp('CV Setting Used: S3 - TARGET');
    end
    disp(' ')

    % run chosen selection method and output CV results
    for ds=[4 3 2 1]
        diary off; diary on;
        disp('--------------------------------------------------------------');
        fprintf('\nData Set: %s\n', datasets{ds});

        % LOAD DATA
        [Y,Sd,St,~,~]=getdata(path,datasets{ds});

        bestaupr = -Inf;
        for sigma=0.25:0.25:2
            for alpha=0:0.1:1
                for eta=0:0.1:1
                    aupr = crossValidation(Y');
                    fprintf('TIME:  %s\t\t', datestr(now));
                    fprintf('sigma=%g\talpha=%g\teta=%g\t\tAUPR:\t%.3g\n', sigma, alpha, eta, aupr)
                    if bestaupr < aupr
                        bestaupr = aupr;
                        bestcomb = [sigma  alpha  eta];
                    end
                end
            end
        end
        disp('Best parameters:')
        disp(['sigma = ' num2str(bestcomb(1))])
        disp(['alpha = ' num2str(bestcomb(2))])
        disp(['  eta = ' num2str(bestcomb(3))])

        disp('--------------------------------------------------------------');
    end
    disp('==============================================================');
end
diary off;

%-------------------------------------------------------------------
%-------------------------------------------------------------------
disp(' ')
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp(' ')
%-------------------------------------------------------------------
%-------------------------------------------------------------------

predictionMethod = 'nbi';

for cvs=1:3
    disp('==============================================================');
    disp(['Prediction method = ' predictionMethod])
    cv_setting = ['S' int2str(cvs)];
    switch cv_setting
        case 'S1', disp('CV Setting Used: S1 - PAIR');
        case 'S2', disp('CV Setting Used: S2 - DRUG');
        case 'S3', disp('CV Setting Used: S3 - TARGET');
    end
    disp(' ')

    % run chosen selection method and output CV results
    for ds=[4 3 2 1]
        diary off; diary on;
        disp('--------------------------------------------------------------');
        fprintf('\nData Set: %s\n', datasets{ds});

        % LOAD DATA
        [Y,Sd,St,~,~]=getdata(path,datasets{ds});

        bestaupr = -Inf;
        for alpha=0:0.1:1
            aupr = crossValidation(Y');
            fprintf('TIME:  %s\t\t', datestr(now));
            fprintf('alpha=%g\t\tAUPR:\t%.3g\n', alpha, aupr)
            if bestaupr < aupr
                bestaupr = aupr;
                bestcomb = alpha;
            end
        end
        disp('Best parameters:')
        disp(['alpha = ' num2str(bestcomb)])

        disp('--------------------------------------------------------------');
    end
    disp('==============================================================');
end
diary off;

%-------------------------------------------------------------------
%-------------------------------------------------------------------
disp(' ')
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp(' ')
%-------------------------------------------------------------------
%-------------------------------------------------------------------

predictionMethod = 'kbmf2k';

for cvs=1:3
    disp('==============================================================');
    disp(['Prediction method = ' predictionMethod])
    cv_setting = ['S' int2str(cvs)];
    switch cv_setting
        case 'S1', disp('CV Setting Used: S1 - PAIR');
        case 'S2', disp('CV Setting Used: S2 - DRUG');
        case 'S3', disp('CV Setting Used: S3 - TARGET');
    end
    disp(' ')

    % run chosen selection method and output CV results
    for ds=[4 3 2 1]
        diary off; diary on;
        disp('--------------------------------------------------------------');
        fprintf('\nData Set: %s\n', datasets{ds});

        % LOAD DATA
        [Y,Sd,St,~,~]=getdata(path,datasets{ds});

        bestaupr = -Inf;
        for rs=10:10:200
            aupr = crossValidation(Y');
            fprintf('TIME:  %s\t\t', datestr(now));
            fprintf('rs=%g\t\tAUPR:\t%.3g\n', rs, aupr)
            if bestaupr < aupr
                bestaupr = aupr;
                bestcomb = rs;
            end
        end
        disp('Best parameters:')
        disp(['rs = ' num2str(bestcomb)])

        disp('--------------------------------------------------------------');
    end
    disp('==============================================================');
end
diary off;

%-------------------------------------------------------------------
%-------------------------------------------------------------------
disp(' ')
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp(' ')
%-------------------------------------------------------------------
%-------------------------------------------------------------------

predictionMethod = 'cmf';
num_iter = 2;
k = 50;

for cvs=1:3
    disp('==============================================================');
    disp(['Prediction method = ' predictionMethod])
    cv_setting = ['S' int2str(cvs)];
    switch cv_setting
        case 'S1', disp('CV Setting Used: S1 - PAIR');
        case 'S2', disp('CV Setting Used: S2 - DRUG');
        case 'S3', disp('CV Setting Used: S3 - TARGET');
    end
    disp(' ')

    % run chosen selection method and output CV results
    for ds=[4 3 2 1]
        diary off; diary on;
        disp('--------------------------------------------------------------');
        fprintf('\nData Set: %s\n', datasets{ds});

        % LOAD DATA
        [Y,Sd,St,~,~]=getdata(path,datasets{ds});

        bestaupr = -Inf;
        for lambda_l=2.^(-3:2)
            for lambda_d=2.^(-4:6)
                for lambda_t=2.^(-4:6)
                    aupr = crossValidation(Y');
                    fprintf('TIME:  %s\t\t', datestr(now));
                    fprintf('lambda_l=%g\tlambda_d=%g\tlambda_t=%g\t\tAUPR:\t%.3g\n', lambda_l, lambda_d, lambda_t, aupr)
                    if bestaupr < aupr
                        bestaupr = aupr;
                        bestcomb = [lambda_l lambda_d lambda_t];
                    end
                end
            end
        end
        disp('Best parameters:')
        disp(['lambda_l = ' num2str(bestcomb(1))])
        disp(['lambda_d = ' num2str(bestcomb(2))])
        disp(['lambda_t = ' num2str(bestcomb(3))])

        disp('--------------------------------------------------------------');
    end
    disp('==============================================================');
end
diary off;

%-------------------------------------------------------------------
%-------------------------------------------------------------------
disp(' ')
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp(' ')
%-------------------------------------------------------------------
%-------------------------------------------------------------------

predictionMethod = 'grmf';
num_iter = 2;
k = 50;

for cvs=1:3
    disp('==============================================================');
    disp(['Prediction method = ' predictionMethod])
    cv_setting = ['S' int2str(cvs)];
    switch cv_setting
        case 'S1', disp('CV Setting Used: S1 - PAIR');
        case 'S2', disp('CV Setting Used: S2 - DRUG');
        case 'S3', disp('CV Setting Used: S3 - TARGET');
    end
    disp(' ')

    % run chosen selection method and output CV results
    for ds=[4 3 2 1]
        diary off; diary on;
        disp('--------------------------------------------------------------');
        fprintf('\nData Set: %s\n', datasets{ds});

        % LOAD DATA
        [Y,Sd,St,~,~]=getdata(path,datasets{ds});

        bestaupr = -Inf;
        for lambda_l=2.^(-5:-1)
            for lambda_d=[0.005 0.01  0.05  0.1  0.15  0.2  0.25  0.3]
                for lambda_t=[0.005 0.01  0.05  0.1  0.15  0.2  0.25  0.3]
                    for p=2:7
                        aupr = crossValidation(Y');
                        fprintf('TIME:  %s\t\t', datestr(now));
                        fprintf('lambda_l=%g\tlambda_d=%g\tlambda_t=%g\tp=%g\t\tAUPR:\t%.3g\n', lambda_l, lambda_d, lambda_t, p, aupr)
                        if bestaupr < aupr
                            bestaupr = aupr;
                            bestcomb = [lambda_l lambda_d lambda_t p];
                        end
                    end
                end
            end
        end
        disp('Best parameters:')
        disp(['lambda_l = ' num2str(bestcomb(1))])
        disp(['lambda_d = ' num2str(bestcomb(2))])
        disp(['lambda_t = ' num2str(bestcomb(3))])
        disp(['       p = ' num2str(bestcomb(4))])

        disp('--------------------------------------------------------------');
    end
    disp('==============================================================');
end
diary off;
