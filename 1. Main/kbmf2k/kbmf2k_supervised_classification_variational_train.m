% Mehmet Gonen (mehmet.gonen@aalto.fi)
% Helsinki Institute for Information Technology HIIT
% Department of Information and Computer Science
% Aalto University School of Science

function state = kbmf2k_supervised_classification_variational_train(Kx, Kz, Y, parameters)
    rand('state', parameters.seed); %#ok<RAND>
    randn('state', parameters.seed); %#ok<RAND>

    Dx = size(Kx, 1);
    Nx = size(Kx, 2);
    Dz = size(Kz, 1);
    Nz = size(Kz, 2);
    R = parameters.R;
    sigmag = parameters.sigmag;

    Lambdax.shape = (parameters.alpha_lambda + 0.5) * ones(Dx, R);
    Lambdax.scale = parameters.beta_lambda * ones(Dx, R);
    Ax.mean = randn(Dx, R);
    Ax.covariance = repmat(eye(Dx, Dx), [1, 1, R]);
    Gx.mean = randn(R, Nx);
    Gx.covariance = eye(R, R);

    Lambdaz.shape = (parameters.alpha_lambda + 0.5) * ones(Dz, R);
    Lambdaz.scale = parameters.beta_lambda * ones(Dz, R);
    Az.mean = randn(Dz, R);
    Az.covariance = repmat(eye(Dz, Dz), [1, 1, R]);
    Gz.mean = randn(R, Nz);
    Gz.covariance = eye(R, R);

    F.mean = (abs(randn(Nx, Nz)) + parameters.margin) .* sign(Y);
    F.covariance = ones(Nx, Nz);    

    KxKx = Kx * Kx';
    KzKz = Kz * Kz';

    lower = -1e40 * ones(Nx, Nz);
    lower(Y > 0) = +parameters.margin;
    upper = +1e40 * ones(Nx, Nz);
    upper(Y < 0) = -parameters.margin;

    lambdax_indices = repmat(logical(eye(Dx, Dx)), [1, 1, R]);
    lambdaz_indices = repmat(logical(eye(Dz, Dz)), [1, 1, R]);

    for iter = 1:parameters.iteration
        if mod(iter, 1) == 0
            fprintf(1, '.');
        end
        if mod(iter, 10) == 0
            fprintf(1, ' %5d\n', iter);
        end

        %%%% update Lambdax
        Lambdax.scale = 1 ./ (1 / parameters.beta_lambda + 0.5 * (Ax.mean.^2 + reshape(Ax.covariance(lambdax_indices), Dx, R)));
        %%%% update Ax
        for s = 1:R
            Ax.covariance(:, :, s) = (diag(Lambdax.shape(:, s) .* Lambdax.scale(:, s)) + KxKx / sigmag^2) \ eye(Dx, Dx);
            Ax.mean(:, s) = Ax.covariance(:, :, s) * (Kx * Gx.mean(s, :)' / sigmag^2);
        end
        %%%% update Gx
        Gx.covariance = (eye(R, R) / sigmag^2 + Gz.mean * Gz.mean' + Nz * Gz.covariance) \ eye(R, R);
        Gx.mean = Gx.covariance * (Ax.mean' * Kx / sigmag^2 + Gz.mean * F.mean');

        %%%% update Lambdaz
        Lambdaz.scale = 1 ./ (1 / parameters.beta_lambda + 0.5 * (Az.mean.^2 + reshape(Az.covariance(lambdaz_indices), Dz, R)));
        %%%% update Az
        for s = 1:R
            Az.covariance(:, :, s) = (diag(Lambdaz.shape(:, s) .* Lambdaz.scale(:, s)) + KzKz / sigmag^2) \ eye(Dz, Dz);
            Az.mean(:, s) = Az.covariance(:, :, s) * (Kz * Gz.mean(s, :)' / sigmag^2);
        end
        %%%% update Gz
        Gz.covariance = (eye(R, R) / sigmag^2 + Gx.mean * Gx.mean' + Nx * Gx.covariance) \ eye(R, R);
        Gz.mean = Gz.covariance * (Az.mean' * Kz / sigmag^2 + Gx.mean * F.mean);
        
        %%%% update F
        output = Gx.mean' * Gz.mean;
        alpha_norm = lower - output;
        beta_norm = upper - output;
        normalization = safenormcdf(beta_norm) - safenormcdf(alpha_norm);
        normalization(normalization == 0) = 1;
        F.mean = output + (safenormpdf(alpha_norm) - safenormpdf(beta_norm)) ./ normalization;
        F.covariance = 1 + (alpha_norm .* safenormpdf(alpha_norm) - beta_norm .* safenormpdf(beta_norm)) ./ normalization - (safenormpdf(alpha_norm) - safenormpdf(beta_norm)).^2 ./ normalization.^2;
    end

    state.Lambdax = Lambdax;
    state.Ax = Ax;
    state.Lambdaz = Lambdaz;
    state.Az = Az;
    state.parameters = parameters;

    state.Ax = rmfield(state.Ax, 'covariance');
    state.Az = rmfield(state.Az, 'covariance');
end