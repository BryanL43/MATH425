% Exercise 1 (Further verification for part d)
fprintf("Exercise 1 (Further verification for part d:\n");
A = rand(4,4);
[Q, R] = qr(A);
disp(Q);
x = rand(4, 1);
disp(norm(Q * x));
disp(norm(x));


% Exercise 2
fprintf("\nExercise 2:\n");
% Part a
nValues = [5]; % Add 10 and 20 later
for n = nValues
    fprintf("For n = %d\n", n);
    H = hilb(n); % Generate Hilbert matrix
    
    disp("--- Part 2A ---")
    [Q, R] = qr(H);
    disp("Q_n:");
    disp(Q);
    disp("R_n:");
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
    x_gauss = H \ b_star;
    disp("Solving Hx = b* using Gaussian elimination:")
    disp(x_gauss);

    % Solving Hx = b* using QR factorization
    [Q, R] = qr(H);
    y = Q' * b_star;
    x_qr = R \ y;
    disp("Solving Hx = b* using QR factorization:")
    disp(x_qr);
end


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

% Following example 4.29 of the textbook to verify
A = [1 1 2; 1 0 -2; -1 2 3];
I = eye(size(A));

%disp("H1:")
v = A(:, 1);
w = I(:, 1);
H1 = myHouseholder(v, w);
%disp(H1);


%disp("H2:")
A2 = H1 * A;
v = A2(:, 2);
v(1) = 0; % Set v_hat2's first element as 0
w = I(:, 2);
H2 = myHouseholder(v, w);
%disp(H2);

%disp("Householder QR:")
Q = H1 * H2;
%disp(Q);
R = H2 * A2;
%disp(R);

% Double check with default MATLAB qr function
%disp("MATLAB QR:")
[q, r] = qr(A);
%disp(q);
%disp(r);
