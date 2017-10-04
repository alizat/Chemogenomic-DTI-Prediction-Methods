% Mehmet Gonen (mehmet.gonen@aalto.fi)
% Helsinki Institute for Information Technology HIIT
% Department of Information and Computer Science
% Aalto University School of Science

function state = kbmf2k_semisupervised_regression_variational_train(Kx, Kz, Y, parameters)
    rand('state', parameters.seed); %#ok<RAND>
    randn('state', parameters.seed); %#ok<RAND>

    Dx = size(Kx, 1);
    Nx = size(Kx, 2);
    Dz = size(Kz, 1);
    Nz = size(Kz, 2);
    R = parameters.R;
    sigmag = parameters.sigmag;
    sigmay = parameters.sigmay;

    Lambdax.shape = (parameters.alpha_lambda + 0.5) * ones(Dx, R);
    Lambdax.scale = parameters.beta_lambda * ones(Dx, R);
    Ax.mean = randn(Dx, R);
    Ax.covariance = repmat(eye(Dx, Dx), [1, 1, R]);
    Gx.mean = randn(R, Nx);
    Gx.covariance = repmat(eye(R, R), [1, 1, Nx]);

    Lambdaz.shape = (parameters.alpha_lambda + 0.5) * ones(Dz, R);
    Lambdaz.scale = parameters.beta_lambda * ones(Dz, R);
    Az.mean = randn(Dz, R);
    Az.covariance = repmat(eye(Dz, Dz), [1, 1, R]);
    Gz.mean = randn(R, Nz);
    Gz.covariance = repmat(eye(R, R), [1, 1, Nz]);

    KxKx = Kx * Kx';
    KzKz = Kz * Kz';

    lambdax_indices = repmat(logical(eye(Dx, Dx)), [1, 1, R]);
    lambdaz_indices = repmat(logical(eye(Dz, Dz)), [1, 1, R]);

    for iter = 1:parameters.iteration
%         if mod(iter, 1) == 0
%             fprintf(1, '.');
%         end
%         if mod(iter, 10) == 0
%             fprintf(1, ' %5d\n', iter);
%         end

        %%%% update Lambdax
        Lambdax.scale = 1 ./ (1 / parameters.beta_lambda + 0.5 * (Ax.mean.^2 + reshape(Ax.covariance(lambdax_indices), Dx, R)));
        %%%% update Ax
        for s = 1:R
            Ax.covariance(:, :, s) = (diag(Lambdax.shape(:, s) .* Lambdax.scale(:, s)) + KxKx / sigmag^2) \ eye(Dx, Dx);
            Ax.mean(:, s) = Ax.covariance(:, :, s) * (Kx * Gx.mean(s, :)' / sigmag^2);
        end
        %%%% update Gx
        for i = 1:Nx
            indices = ~isnan(Y(i, :));
            Gx.covariance(:, :, i) = (eye(R, R) / sigmag^2 + (Gz.mean(:, indices) * Gz.mean(:, indices)' + sum(Gz.covariance(:, :, indices), 3)) / sigmay^2) \ eye(R, R);
            Gx.mean(:, i) = Gx.covariance(:, :, i) * (Ax.mean' * Kx(:, i) / sigmag^2 + Gz.mean(:, indices) * Y(i, indices)' / sigmay^2);
        end

        %%%% update Lambdaz
        Lambdaz.scale = 1 ./ (1 / parameters.beta_lambda + 0.5 * (Az.mean.^2 + reshape(Az.covariance(lambdaz_indices), Dz, R)));
        %%%% update Az
        for s = 1:R
            Az.covariance(:, :, s) = (diag(Lambdaz.shape(:, s) .* Lambdaz.scale(:, s)) + KzKz / sigmag^2) \ eye(Dz, Dz);
            Az.mean(:, s) = Az.covariance(:, :, s) * (Kz * Gz.mean(s, :)' / sigmag^2);
        end
        %%%% update Gz
        for j = 1:Nz
            indices = ~isnan(Y(:, j));
            Gz.covariance(:, :, j) = (eye(R, R) / sigmag^2 + (Gx.mean(:, indices) * Gx.mean(:, indices)' + sum(Gx.covariance(:, :, indices), 3)) / sigmay^2) \ eye(R, R);
            Gz.mean(:, j) = Gz.covariance(:, :, j) * (Az.mean' * Kz(:, j) / sigmag^2 + Gx.mean(:, indices) * Y(indices, j) / sigmay^2);
        end
    end

    state.Lambdax = Lambdax;
    state.Ax = Ax;
    state.Lambdaz = Lambdaz;
    state.Az = Az;
    state.parameters = parameters;

    state.Ax = rmfield(state.Ax, 'covariance');
    state.Az = rmfield(state.Az, 'covariance');
end