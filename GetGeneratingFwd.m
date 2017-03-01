function [Ggen, XYZGenAct] = GetGeneratingFwd(XYZGen, G2d, R)
    % Assign topographies to generator coordinates(do not recompute, just find the closest one)
    for i=1:size(XYZGen, 1)
        d = repmat(XYZGen(i,:), size(R,1), 1) - R;
        d = sum(d .* d, 2);
        [~, ind] = min(d);
        XYZGenAct(i,:) = R(ind,:);
        Ggen(:, i) = G2d(:, ind * 2 - 1); % take the first dipole in the tangent plane
        GenInd(i) = ind; 
    end;
end
