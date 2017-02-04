function aucs = calc_auc(X,Y)
    [n_monte, n_snr, n_pts] = size(X);
    aucs = zeros(n_monte, n_snr);
    for i = 1:n_monte
        for j = 1:n_snr
            aucs(i,j) = trapz(squeeze(X(i,j,:)), squeeze(Y(i,j,:)));
        end
    end
end

