% Exercise 1
disp("Exercise 1:");
A = [1 0 1; 2 1 0; -1 -2 3; 0 -1 2];
b = [1; 1; 2; -2];

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
x_qrstar = myBackwardSubstitution(R, Q' * b);
disp("Least square using QR-factorization:");
disp(x_qrstar);

% Exercise 3
disp("Exercise 3:");
years = 1989:1999;
prices = [86.4, 89.8, 92.8, 96.0, 99.6, 103.1, 106.3, 109.5, 113.3, 120.0, 129.5];
x = years - 1989; % Use base to turn years into 0, 1, ..., 10

% Construct A matrix where 1st col = constant 1; 2nd col = years starting
% from 0
A = [ones(length(x), 1) x'];
b = prices'; % transposed to vertical vector

x_star = (A' * A)^-1 * A' * b;

% Acquire equation from x_star
alpha = x_star(1);
beta = x_star(2);
fprintf("Equation of least square line is about: y = %.2f + %.2fx\n", alpha, beta);

% Estimate the median price of a house in the year 2005 & 2010
x_2005 = 2005 - 1989;
x_2010 = 2010 - 1989;
price_2005 = alpha + x_2005 * beta;
price_2010 = alpha + x_2010 * beta;

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

% Part b:
zeta_8 = exp(1i * 2 * pi / n); % ζ8
omega = zeros(n, n); % matrix to hold w_ks

for k = 0:n-1
    for j = 0:n-1
        omega(j + 1, k + 1) = zeta_8^(j * k); % ζ8
    end
end

% Part c:
c = zeros(n, 1); % init coefficient vector

for k = 1:n
    c(k) = (1 / n) * dot(omega(:, k), f);  % Compute the inner product and scale by 1/8
end

% Part d:
x = linspace(0, 2*pi, 100); % range of values from 0 to 2*pi

%p1 = real(c(1)) * ones(size(x)); % % real part init with 1st term
p1 = zeros(size(x));
for k = 0:n-1
    p1 = p1 + real(c(k + 1)) * cos(k * x) - imag(c(k + 1)) * sin(k * x);
end

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

% Labels
title("Comparison of f(x) = x^2 and p_1(x) via Fourier approximation");
xlabel("x");
ylabel("f(x)");
legend("f(x) = x^2", "p_1(x)");

% Part f:
zeta_8 = exp(1i * 2 * pi / n); % ζ8
omega = zeros(n, n); % matrix to hold w_ks

% Iterate over k from -4 to 3
i = 1;
for k = -4:3
    for j = 0:n-1
        omega(j + 1, i) = zeta_8^(j * k); % ζ8
    end
    i = i + 1;
end

c = zeros(n, 1); % init modified coefficient vector

for k = 1:n
    c(k) = (1 / n) * dot(omega(:, k), f);  % Compute the inner product and scale by 1/8
end

% Part g:
q1 = zeros(size(x));

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

% Labels
title("Comparison of f(x) = x^2 and q_1(x) via Fourier reconstruction");
xlabel("x");
ylabel("f(x)");
legend("f(x) = x^2", "q_1(x)");
