% Imported functions to support Exercise 2
function U = myGaussianElimination(A)
    % Check if matrix is an n*(n+1) matrix
    [n, m] = size(A);
    if m ~= n+1
        error("Input matrix must be an n*(n+1) matrix.");
    end

    for j = 1: n % Iterate through rows
        % Check for zero pivot
        if A(j, j) == 0
            error("Pivot element is zero.");
        end

        for i = j + 1: n % Iterate through column
            factor = -1 * (A(i, j) / A(j, j)); % Compute elimiantion factor

            % Apply row operation to make the entries below the pivot zero
            A(i, j:end) = A(i, j:end) + factor * A(j, j:end);
        end
    end
    U = A;
end

function x = myBackwardSubstitution(U, c)
    [n_U, m_U] = size(U);
    [n_c, m_c] = size(c);

    % Check matrices dimensions
    if n_c ~= n_U || m_U ~= n_U || m_c ~= 1
        error("Dimensions of U is not n*n or c is not n*1.");
    end

    x = zeros(n_U, 1); % Initialize solution vector
    
    % Iterate through the rows in ascending order
    for i = n_U: -1: 1
        sum = c(i); % Acquire RHS values for the current row
        
        % Subtract the already computed x values
        for j = i + 1: n_U
            sum = sum - U(i, j) * x(j);
        end

        % Divide the diagonal element
        x(i) = sum / U(i, i);
    end
end

function solution = myLinearSolution(A, b)
    [n_A, m_A] = size(A);
    [n_b, m_b] = size(b);
    
    % Check matrices dimensions
    if n_b ~= n_A || m_A ~= n_A || m_b ~= 1
        error("Dimensions of A is not n*n or b is not n*1.");
    end

    Ab = [A b]; % Augment A & b
    
    upperTriMatrix = myGaussianElimination(Ab);
    
    % Split the upper-triangular matrix from the augmented RHS
    U = upperTriMatrix(:, 1: end - 1);
    c = upperTriMatrix(:, end);

    solution = myBackwardSubstitution(U, c);
end

% Exercise 2
fprintf("\nExercise 2:\n");
% Part a
nValues = [5, 10, 20];
for n = nValues
    fprintf("For n = %d\n", n);
    H = hilb(n); % Generate Hilbert matrix
    
    disp("--- Part 2A ---")
    [Q, R] = qr(H);
    fprintf("Q_%d\n", n);
    disp(Q);
    fprintf("R_%d\n", n);
    disp(R);
    
    disp("--- Part 2B ---")
    % Define x*
    x_star = zeros(n, 1);
    for i = 1:n
        x_star(i) = (-1)^i * (i / (i + 1));
    end
    disp("x*:");
    disp(x_star);

    % Compute b* = Hx*
    b_star = H * x_star;
    disp("b*:");
    disp(b_star);

    % Solving Hx = b* using first Gaussian elimination
    x_gauss_star = myLinearSolution(H, b_star);
    disp("Solving Hx = b* using Gaussian elimination:");
    disp("x_gauss* = ");
    disp(x_gauss_star);

    % Solving Hx = b* using QR factorization
    [Q, R] = qr(H);
    y = Q' * b_star; % Computing LHS: Rx = Q^T * b*
    x_qr_star = myBackwardSubstitution(R, y); % Solving Rx = y
    disp("Solving Hx = b* using QR factorization:");
    disp("x_qr* = ");
    disp(x_qr_star);

    disp("Error margins:");
    disp("x* - x_gauss* = ");
    disp(norm(x_star - x_gauss_star));

    disp("x* - x_qr* = ");
    disp(norm(x_star - x_qr_star));
end

disp("--- Part 2C ---")
fprintf("Comparing the two methods to the correct solution x*:\n" + ...
    "QR-factorization is more stable then Gaussian elimination as the Hilbert matrix gets larger, shown by the error magins\n" + ...
    "calculated by the difference of the correct solution with the result of QR/Gaussian.\n\n" + ...
    "The pros of using Gaussian Elimination is that it is less complex, requires less computation, and easy to implement.\n" + ...
    "It is also efficient for small matrices. However, Gaussian is more numerically instable especially for Hilbert matrices.\n" + ...
    "It also fails if a pivot is zero.\n\n" + ...
    "The pros of QR-factorization is that is is more numerically stable especially for large and Hilbert matrices.\n" + ...
    "However, it requires more complex and requires additional computation compared to Gaussian elimination.\n")

% Exercise 3
fprintf("\nExercise 3:\n");
function H = myHouseholder(v, w)
    % Ensure non-zero vectors
    if norm(v) == 0 || norm(w) == 0
        error("v or w is a zero vector.");
    end
    
    % Ensure v and w has the same dimension in R^n
    if length(v) ~= length(w)
        error("v and w is not both in R^n.");
    end
    
    % Normalize the vectors v and w
    v_hat = v / norm(v);
    w_hat = w / norm(w);
    
    % Compute the Householder matrix
    u = v_hat - w_hat;
    I = eye(length(v));
    H = I - 2 * (u * u') / (norm(u)^2);
end

% Testing myHouseholder function with 3 random vectors in R^4
v1 = randn(4, 1)
v2 = randn(4, 1)
v3 = randn(4, 1);

disp("Test 1 (myHouseholder(v1, v2)):");
v1_hat = v1 / norm(v1);
v2_hat = v2 / norm(v2);
H1 = myHouseholder(v1, v2);

disp('H1 * v1_hat and v2_hat:');
disp([H1 * v1_hat, v2_hat]);

disp("Test 2 (myHouseholder(v2, v3)):");
v3_hat = v3 / norm(v3);
H2 = myHouseholder(v2, v3);

disp('H2 * v2_hat and v3_hat:');
disp([H2 * v2_hat, v3_hat]);

disp("Test 3 (myHouseholder(v1, v3)):");
H3 = myHouseholder(v1, v3);

disp('H3 * v1_hat and v3_hat:');
disp([H3 * v1_hat, v3_hat]);


% Following example 4.29 of the textbook to verify integrity
disp("Verifying myHouseholder function integrity with textbook example 4.29:")
A = [1 1 2; 1 0 -2; -1 2 3];
I = eye(size(A));

disp("H1:")
v = A(:, 1);
w = norm(A(:, 1)) * I(:, 1);
H1 = myHouseholder(v, w);
disp(H1);

disp("H2:")
A2 = H1 * A;
v = A2(:, 2);
v(1) = 0; % Set v_hat2's first element as 0
w = norm(v) * I(:, 2);
H2 = myHouseholder(v, w);
disp(H2);

disp("Householder QR:")
Q = H1 * H2;
disp(Q);
R = H2 * A2;
disp(R);
