function A = make_A(N, alpha, h)
    % Input: N - number of points
    %        alpha - constant in the equation
    %        h - grid spacing
    % Output: A - NxN matrix

    % Create a diagonal matrix with -2(alpha)-2/h^2 as the main diagonal
    main_diag = (-2 * alpha - 2 / h^2) * ones(N, 1);
    A = diag(main_diag);

    % Create the off-diagonal entries with 1/h^2
    off_diag = 1 / h^2 * ones(N - 1, 1);

    % Add the off-diagonal entries to A
    A = A + diag(off_diag, 1) + diag(off_diag, -1);

    % Set the top-right and bottom-left corner entries for periodic boundary conditions
    A(1, N) = 1 / h^2;
    A(N, 1) = 1 / h^2;
end