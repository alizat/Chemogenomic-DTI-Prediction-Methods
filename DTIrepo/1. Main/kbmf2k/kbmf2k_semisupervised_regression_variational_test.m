% Mehmet Gonen (mehmet.gonen@aalto.fi)
% Helsinki Institute for Information Technology HIIT
% Department of Information and Computer Science
% Aalto University School of Science

function prediction = kbmf2k_semisupervised_regression_variational_test(Kx, Kz, state)
    Nx = size(Kx, 2);
    Nz = size(Kz, 2);
    R = size(state.Ax.mean, 2);

    prediction.Gx.mean = zeros(R, Nx);
    for s = 1:R
        prediction.Gx.mean(s, :) = state.Ax.mean(:, s)' * Kx;
    end

    prediction.Gz.mean = zeros(R, Nz);
    for s = 1:R
        prediction.Gz.mean(s, :) = state.Az.mean(:, s)' * Kz;
    end

    prediction.Y.mean = prediction.Gx.mean' * prediction.Gz.mean;
end