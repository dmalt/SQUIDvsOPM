function min_dist = ExamineGrid(R)
    n_src = size(R, 1);

    for i = 1:n_src
        temp_loc = sum( ( R - repmat(R(i,:), [size(R, 1), 1]) ) .^ 2, 2  );
        min_dist(i) = sqrt(min(temp_loc(temp_loc > 0)));
    end
    figure;
    histogram(min_dist, 50);
end
