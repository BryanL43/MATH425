% Exercise #1
A = [2 -1 2; -1 -1 3; 3 0 -2];
b = [2; 1; 1];

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

myLinearSolution(A, b)


% Exercise #2
A = [-8 -2 3 1; 1 -2 0 2; -4 -1 3 2; 4 1 -1 -1];

% Acquire LU factorization of A
[L, U, I] = lu(A)

% Part h
b = [-2; 6; -5; 1];

% 1st step: solve for y in Ly = b using forward substitution
y = L \ (I * b);

% 2nd step: solve for x in Ux = y
x = U \ y;

solu = myLinearSolution(A, b);

disp("Is solving x using LU factorization of A equals to using myLinearSolution (1 = True; 0 = false): ");
disp(isequal(x, solu));


% Exercise #4
H_5 = hilb(5)
H_10 = hilb(10);
H_15 = hilb(15);
H_20 = hilb(20);

Inverse_H_5 = inv(H_5)
Inverse_H_10 = inv(H_10)
Inverse_H_15 = inv(H_15)
Inverse_H_20 = inv(H_20)

% The feature that jumps at me is the large number (the ...e+... numbers) multiplied with the
% matrix. It even starts displaying a warning of inaccuracy at hilb(15).

x = rand(15, 1)
b = H_15 * x

x_computed = Inverse_H_15 * b

check = H_15 * x_computed

% Using x_computed is quite close to just x but with a margin of error.
% There also seems to be some negative values in the x_computed compared to
% regular x.
