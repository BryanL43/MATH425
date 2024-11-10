% Exercise 1
disp("Exercise 1:");
b = [1; 1; 2; -2];

v_1 = [1; 2; -1; 0];
v_2 = [0; 1; -2; -1];
v_3 = [1; 0; 3; 2];
A = [v_1 v_2 v_3];

% Augment A' * A and A' * b to solve via row operation
A_hat = A' * A;
b_hat = A' * b;
A_tilde = [A_hat b_hat];
disp(rref(A_tilde)); % Note: A is not full-rank. so let z be free variable

% The computed solution vector
z = 0; % arbitrary value due to linear dependence
y = -0.5 + 2 * z;
x = 0.5 - z;
v = [x; y; z];
disp(v);

% Now, computing the closest points from b vector to the subspace spanned
% by v_1, v_2, and v_3
closestPoint = A * v;
fprintf("The closest points is:\n");
disp(closestPoint); % vertical vector of all 0.5

% Answer integrity check via direct least square:
fprintf("Answer check via least square:\n");
x_star = A \ b;
closestPoint = A * x_star;
disp(closestPoint);

% Exercise 2
disp("Exercise 2:");
A = [1 2 -1; 0 -2 3; 1 5 -1; -3 1 1];
b = [0; 5; 6; 8];

disp("Calculate least square via x* = (A^T * A)^-1 * A^T * b:");
x_star = (A' * A)^-1 * A' * b;
disp(x_star);

% Imported backward-sub function from previous assignments
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

[Q, R] = qr(A, 0); % "Economy-sized" QR for optimized use on A with more rows than cols
x_qrstar = myBackwardSubstitution(R, Q' * b); % Solving Rx = Q^T * b
disp("Least square using QR-factorization:");
disp(x_qrstar);

% Exercise 3
disp("Exercise 3:");
years = 1989:1999;
prices = [86.4, 89.8, 92.8, 96.0, 99.6, 103.1, 106.3, 109.5, 113.3, 120.0, 129.5];

% Construct A matrix where 1st col = constant 1; 2nd col = years starting
% from 0
A = [ones(length(years), 1) years'];
b = prices'; % transposed to vertical vector

x_star = (A' * A)^-1 * A' * b;

% Acquire equation from x_star
alpha = x_star(1);
beta = x_star(2);
fprintf("Equation of least square line is about: y = %.2f + %.2fx\n", alpha, beta);

% Estimate the median price of a house in the year 2005 & 2010
price_2005 = alpha + 2005 * beta;
price_2010 = alpha + 2010 * beta;

fprintf("Estimated median price in 2005: $%.2f thousand\n", price_2005);
fprintf("Estimated median price in 2010: $%.2f thousand\n\n", price_2010);

% Exercise 4
disp("Exercise 4:");

% Part a:
n = 8; % # of sample points
f = zeros(n, 1); % init sample vector f

% Compute f(j * (2pi)/8) for j = 0, 1, ..., 7
for j = 0:n-1
    x = (j * 2 * pi) / n;
    f(j + 1) = x^2; % solves f(x) & accounts for matlab indexing
end
disp("Sample vector f:");
disp(f);

% Part b:
zeta_8 = exp(1i * 2 * pi / n); % ζ8
omega = zeros(n, n); % matrix to hold w_ks

for k = 0:n-1
    for j = 0:n-1
        omega(j + 1, k + 1) = zeta_8^(j * k); % ζ8^kj
    end
end
disp("omega vectors:");
disp(omega);

% Part c:
c = zeros(n, 1); % init coefficient vector

for k = 1:n
    % Compute c_k = ⟨f,ω_k⟩ for k = 0,...,7 and scale by 1/8
    c(k) = (1 / n) * dot(omega(:, k), f);
end
disp("c_k vector:");
disp(c);

% Part d:
% Get range of values from 0 to 2*pi to smooth the graph out
x = linspace(0, 2*pi, 100);

% Compute p1(x) from p(x) = p1(x) + ip2(x).
% Let a be the real coefficient of c and b be the imaginary coefficent of c.
% p(x) is computed as p(x) = (a + ib)(cos(nx) + isin(nx)), which condenses
% into p(x) = (acos(nx) - bsin(nx)) + i(asin(nx) + bcos(nx)).
% p1(x) = acos(nx) - bsin(nx), p2(x) = aisin(nx) + bisin(nx).
p1 = zeros(size(x));
for k = 0:n-1
    p1 = p1 + real(c(k + 1)) * cos(k * x) - imag(c(k + 1)) * sin(k * x);
end
disp("p_1 values:");
disp(p1);

% Part e:
f_x = x.^2; % f(x) = x^2

% Plot f(x) and p_1(x)
figure;
plot(x, f_x, "r-", "LineWidth", 1.5);
hold on;
plot(x, p1, 'b--', 'LineWidth', 1.5);
hold off;

xticks([0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
grid on;

% Labels
title("Comparison of f(x) = x^2 and p_1(x) via Fourier approximation");
xlabel("x");
ylabel("f(x)");
legend("f(x) = x^2", "p_1(x)");

% Part f:
zeta_8 = exp(1i * 2 * pi / n); % ζ8
omega = zeros(n, n); % matrix to hold w_ks

% Iterate over k from -4 to 3
i = 1; % for index mapping
for k = -4:3
    for j = 0:n-1
        omega(j + 1, i) = zeta_8^(j * k); % ζ8^(kj)
    end
    i = i + 1;
end
disp("Shifted omega vectors:");
disp(omega);

c = zeros(n, 1); % init modified coefficient vector

for k = 1:n
    c(k) = (1 / n) * dot(omega(:, k), f);  % Compute the inner product and scale by 1/8
end
disp("Shifted c_k vectors:");
disp(c);

% Part g:
q1 = zeros(size(x));

% Same explanation for computing q1 as explained in part d. We only want to
% graph the real numbers.
i = 1; % To keep track of c vector position to read from
for k =-4:3
    q1 = q1 + real(c(i)) * cos(k * x) - imag(c(i)) * sin(k * x);
    i = i + 1;
end

% Part h:
figure;
plot(x, f_x, "r-", "LineWidth", 1.5);
hold on;
plot(x, q1, 'b--', 'LineWidth', 1.5);
hold off;

xticks([0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
grid on;

% Labels
title("Comparison of f(x) = x^2 and q_1(x) via Fourier reconstruction");
xlabel("x");
ylabel("f(x)");
legend("f(x) = x^2", "q_1(x)");
