function [Cs, Ps, IND] = iDICS_1D(C, G)

    Nch = size(C, 1);
    iC = inv(C + 10 * trace(C) / Nch * eye(Nch));
    Ns = fix(size(G, 2)); % assume tangent space dimension of 2

    A = zeros(size(G'));
    for i=1:Ns
        L = G(:, i);
        A(i, :) = inv(L' * iC * L) * L' * iC;
    end

    ACA = A * C * A';
    Ps = diag(ACA);
    dACA_recip = 1. / sqrt(diag((ACA)));

    % for i=1:size(ACA,1)
    %     ACA(:,i) = ACA(:,i).*(dACA_recip);
    % end
    % 
    % for i=1:size(ACA,1)
    %     ACA(i,:) = ACA(i,:).*(dACA_recip');
    % end

    %extract the upper triangle(exclude the diagonal)
    % and make sure it sorted naturally(1st row, second row, etc)
    T = triu(ones(size(ACA)), 1);
    Cs = abs(imag(ACA(T == 1)));
    [I, J] = ind2sub(size(ACA), find(T == 1));
    [~, key] =  sort(I, 'descend');
    IND = [I(key), J(key)];
    Cs = Cs(key);

    fprintf('\n Done\n');
end
