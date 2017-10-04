function getParameters(classifier,cv_setting,ds)
%
% This function retrieves optimal parameter values based on the supplied
% algorithm, the CV scenario used and the dataset on which prediction will
% be performed.
%
% INPUT:
%  classifier:  algorithm to be used for DTI prediction
%  cv_setting:  cross validation setting used (S1/S2/S3)
%  ds:          dataset (4:NR, 3:GPCR, 2:IC, 1:E)
%
% OUTPUT:
%  params:   optimal parameter values for given inputs
%

    global num_iter k lambda_l lambda_d lambda_t p sigma alpha eta rs regcoef lambda

    switch classifier

        % -----------------------------------------

        case 'grmf'
            num_iter = 2;
            k = 100;
            switch cv_setting

                case 'S1'
                    switch(ds)
                        case 4
                            lambda_l = 0.125;
                            lambda_d = 0.15;
                            lambda_t = 0.05;
                                   p = 2;

                        case 3
                            lambda_l = 0.5;
                            lambda_d = 0.3;
                            lambda_t = 0.3;
                                   p = 3;

                        case 2
                            lambda_l = 0.5;
                            lambda_d = 0.3;
                            lambda_t = 0.3;
                                   p = 3;

                        case 1
                            lambda_l = 0.5;
                            lambda_d = 0.3;
                            lambda_t = 0.3;
                                   p = 7;
                    end

                case 'S2'
                    switch(ds)
                        case 4
                            lambda_l = 0.0625;
                            lambda_d = 0.05;
                            lambda_t = 0.1;
                            p = 2;

                        case 3
                            lambda_l = 0.0625;
                            lambda_d = 0.05;
                            lambda_t = 0.01;
                            p = 7;

                        case 2
                            lambda_l = 0.25;
                            lambda_d = 0.15;
                            lambda_t = 0.01;
                            p = 6;

                        case 1
                            lambda_l = 0.25;
                            lambda_d = 0.2;
                            lambda_t = 0.2;
                            p = 4;
                    end
                    
                case 'S3'
                    switch(ds)
                        case 4
                            lambda_l = 0.125;
                            lambda_d = 0.2;
                            lambda_t = 0.05;
                            p = 2;

                            %best result for NR dataset under S3:
                            % k = 100;
                            % lambda_l = 0.125;  % SAME
                            % lambda_d = 0.5;
                            % lambda_t = 0.05;   % SAME
                            % p_d = 8;
                            % p_t = 1;

                        case 3
                            lambda_l = 0.0625;
                            lambda_d = 0.3;
                            lambda_t = 0.01;
                            p = 2;

                        case 2
                            lambda_l = 0.25;
                            lambda_d = 0.2;
                            lambda_t = 0.2;
                            p = 4;

                        case 1
                            lambda_l = 0.25;
                            lambda_d = 0.2;
                            lambda_t = 0.2;
                            p = 5;

                    end
                    
            end
            fprintf('lambda_l=%g, lambda_d=%g, lambda_t=%g, p=%g\n', lambda_l, lambda_d, lambda_t, p)

        % -----------------------------------------
            
        case 'cmf'
            num_iter = 2;
            k = 100;
            switch cv_setting

                case 'S1'
                    switch(ds)
                        case 4
                            lambda_l = 1;
                            lambda_d = 2;
                            lambda_t = 0.125;

                        case 3
                            lambda_l = 1;
                            lambda_d = 0.25;
                            lambda_t = 0.5;

                        case 2
                            lambda_l = 1;
                            lambda_d = 0.25;
                            lambda_t = 0.5;

                        case 1
                            lambda_l = 1;
                            lambda_d = 0.0625;
                            lambda_t = 2;
                    end
                    
                case 'S2'
                    switch(ds)
                        case 4
                            lambda_l = 1;
                            lambda_d = 4;
                            lambda_t = 32;

                        case 3
                            lambda_l = 2;
                            lambda_d = 8;
                            lambda_t = 0.125;

                        case 2
                            lambda_l = 4;
                            lambda_d = 32;
                            lambda_t = 0.25;

                        case 1
                            lambda_l = 2;
                            lambda_d = 64;
                            lambda_t = 0.125;
                    end
                    
                case 'S3'
                    switch(ds)
                        case 4
                            lambda_l = 0.25;
                            lambda_d = 0.25;
                            lambda_t = 32;

                        case 3
                            lambda_l = 0.5;
                            lambda_d = 0.0625;
                            lambda_t = 64;

                        case 2
                            lambda_l = 0.5;
                            lambda_d = 0.0625;
                            lambda_t = 32;

                        case 1
                            lambda_l = 1;
                            lambda_d = 0.0625;
                            lambda_t = 2;
                    end
                    
            end
            fprintf('lambda_l=%g, lambda_d=%g, lambda_t=%g\n', lambda_l, lambda_d, lambda_t)

        % -----------------------------------------
            
        case 'rls_wnn'
            switch cv_setting

                case 'S1'
                    switch(ds)
                        case 4
                            sigma = 0.5;
                            alpha = 0.2;
                              eta = 0.6;

                        case 3
                            sigma = 1;
                            alpha = 0.5;
                              eta = 0.2;

                        case 2
                            sigma = 2;
                            alpha = 0.4;
                              eta = 0.1;

                        case 1
                            sigma = 2;
                            alpha = 0.5;
                              eta = 0.1;
                    end

                case 'S2'
                    switch(ds)
                        case 4
                            sigma = 0.5;
                            alpha = 0.5;
                            eta = 0.6;

                        case 3
                            sigma = 0.5;
                            alpha = 1;
                            eta = 0.7;

                        case 2
                            sigma = 0.25;
                            alpha = 1;
                            eta = 0.7;

                        case 1
                            sigma = 0.5;
                            alpha = 1;
                            eta = 0.7;
                    end
                    
                case 'S3'
                    switch(ds)
                        case 4
                            sigma = 0.25;
                            alpha = 0.1;
                            eta = 0.4;

                        case 3
                            sigma = 0.5;
                            alpha = 0.9;
                            eta = 0.4;

                        case 2
                            sigma = 2;
                            alpha = 0.6;
                            eta = 0.1;

                        case 1
                            sigma = 0.5;
                            alpha = 1;
                            eta = 0.6;
                    end
                    
            end
            fprintf('sigma=%g, alpha=%g, eta=%g\n', sigma, alpha, eta)

        % -----------------------------------------
            
        case 'nbi'
            alpha = 0.5;
            fprintf('alpha=%g\n', alpha)

        % -----------------------------------------

        case 'kbmf2k'
            switch cv_setting

                case 'S1'
                    switch(ds)
                        case 4
                            rs = 30;

                        case 3
                            rs = 100;

                        case 2
                            rs = 80;

                        case 1
                            rs = 170;
                    end

                case 'S2'
                    switch(ds)
                        case 4
                            rs = 50;

                        case 3
                            rs = 100;   %170

                        case 2
                            rs = 100;

                        case 1
                            rs = 200;           % aupr = 0.2542 (0.00976132)
                    end
                    
                case 'S3'
                    switch(ds)
                        case 4
                            rs = 120;   %20

                        case 3
                            rs = 190;           % aupr = 0.527728

                        case 2
                            rs = 130;   %170

                        case 1
                            rs = 200;           % aupr = 0.671863 (0.0240449)
                    end
                    
            end
            fprintf('rs=%g\n', rs)

        % -----------------------------------------
            
        case 'kron_rls_mkl'
            lambda = 1;
            regcoef = 0.25;

    end

end