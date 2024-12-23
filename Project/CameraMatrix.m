% @brief: Computes the similarity transformation matrix T with its
% respective normalized image points and the similarity transformation
% matrix U with its respective normalized world points.
% 
% @param x: 2D image points, a matrix of size (n x 3) where each row is
% [x, y, 1]
% @param X: 3D world points, a matrix of size (n x 4) where each row is
% [X, Y, Z, 1]
% return [T, x_norm, U, X_norm]: The two similarity trasnformation matrices
% and the normalized image and world points.
function [T, x_norm, U, X_norm] = normalizePoints(x, X)
    % Ensure there are n >= 6 points to compute camera matrix P
    X_n = size(X, 1);
    x_n = size(x, 1);
    if X_n < 6 || x_n < 6
        error("At least 6 points (both X & x) are required to compute the camera matrix.");
    end

    % Ensure 3D world points matrix has 4 columns &
    % 2D image points matrix has 3 columns
    if size(X, 2) ~= 4 || size(x, 2) ~= 3
        error("Incorrect dimensions. X should have 4 columns & x should have 3 columns.");
    end
    
    % Acquire the center point of all mapped x (2D) points (centroid)
    centroid_x = mean(x(:, 1:2));
    
    % Set centroid_x's coordinate origin to (0, 0)^T
    x_centered = x(:, 1:2) - centroid_x;
    
    % Compute x_centered's average distance from origin
    avg_dist_x = mean(sqrt(sum(x_centered.^2, 2)));

    % Scale the points to make the average distance from origin to sqrt(2)
    scale_x = sqrt(2) / avg_dist_x;

    % Create the image points transformation matrix T
    T = [scale_x, 0, -scale_x * centroid_x(1);
         0, scale_x, -scale_x * centroid_x(2);
         0, 0, 1];

    % Normalize 2D points with T
    x_norm = (T * x')';

    % Acquire the center point of all mapped X (3D) points (centroid)
    centroid_X = mean(X(:, 1:3));

    % Set centroid_X's coordinate origin to (0, 0)^T
    X_centered = X(:, 1:3) - centroid_X;
    
    % Compute X_centered's average distance from origin
    avg_dist_X = mean(sqrt(sum(X_centered.^2, 2)));
    
    % Scale the points to make the average distance from origin to sqrt(2)
    scale_X = sqrt(3) / avg_dist_X;
    
    % Create the space points transformation matrix U
    U = [scale_X, 0, 0, -scale_X * centroid_X(1);
         0, scale_X, 0, -scale_X * centroid_X(2);
         0, 0, scale_X, -scale_X * centroid_X(3);
         0, 0, 0, 1];

    % Normalize 3D points with U
    X_norm = (U * X')';
end

% @brief: Computes the camera matrix P given a set of 6 or more 3D points
% with its respective 2D point on an image plane for an uncalibrated
% camera. Assuming an affine camera model non-respective to perspective.
% 
% @param x: 2D image points, a matrix of size (n x 3) where each row is
% [x, y, 1]
% @param X: 3D world points, a matrix of size (n x 4) where each row is
% [X, Y, Z, 1]
% return P: The computed camera projection matrix (3 x 4)
function P = computeCameraMatrix(x, X)
    % Ensure there are n >= 6 points to compute camera matrix P
    X_n = size(X, 1);
    x_n = size(x, 1);
    if X_n < 6 || x_n < 6
        error("At least 6 points (both X & x) are required to compute the camera matrix.");
    end

    % Ensure 3D world points matrix has 4 columns &
    % 2D image points matrix has 3 columns
    if size(X, 2) ~= 4 || size(x, 2) ~= 3
        error("Incorrect dimensions. X should have 4 columns & x should have 3 columns.");
    end

    n = X_n;

    % Construct the matrix A
    A = zeros(2 * n, 12); % Each point has corresponding (a_xi)^T & (a_yi)^T

    for i = 1:n
        X_i = X(i, :); % Current 3D point [X, Y, Z, 1] (X_i)
        u_i = x(i, 1); % Image point x-coordinate (x_i)
        v_i = x(i, 2); % Image point y-coordinate (y_i)
        
        % Compute the corresponding (a_xi)^T row:
        % (a_xi)^T = [-X_i, -Y_i, -Z_i, -1, 0, 0, 0, 0, (y_i * X_i), (y_i *
        % Y_i), (y_i * Z_i), y_i]
        A(2 * i - 1, :) = [-X_i, zeros(1, 4), u_i * X_i];
        
        % Compute the corresponding (a_yi)^T row:
        % (a_xi)^T = [0, 0, 0, 0, -X_i, -Y_i, -Z_i, -1, (x_i * X_i), (x_i *
        % Y_i), (x_i * Z_i), x_i]
        A(2 * i, :) = [zeros(1, 4), -X_i, v_i * X_i];
    end
    
    % Solving for Ap = 0 using Singular Value Decomposition (SVD)
    [~, ~, Q] = svd(A);
 
    % Take the last column (smallest singular vector) of Q as solution to minimize σ = Σ_ii
    p = Q(:, end);
 
    % Reshape p into the 3x4 camera matrix P
    P = reshape(p, 4, 3)';
end

% Define 3D world points (m x 4), homogeneous coordinates
X = [0, 0, 0, 1;
     1, 0, 0, 1;
     1, 1, 0, 1;
     0, 1, 0, 1;
     0, 0, 1, 1;
     1, 0, 1, 1];

% Define corresponding 2D image points (m x 3), homogeneous coordinates
x = [100, 200, 1;
     200, 200, 1;
     200, 300, 1;
     100, 300, 1;
     120, 220, 1;
     220, 220, 1];

[T, x_normalized, U, X_normalized] = normalizePoints(x, X);

P_norm = computeCameraMatrix(x_normalized, X_normalized);
disp("Normalized Camera Projection Matrix P:")
disp(P_norm);

% Acquire original camera matrix by denormalize the normalized camera matrix
P = inv(T) * P_norm * U;
disp("Final camera matrix P:");
disp(P);

% Verification: Project the 3D world point's 2D image point using the final camera matrix
x_projected = (P * X')';
for i = 1:size(x_projected, 1)
    x_projected(i, :) = x_projected(i, :) / x_projected(i, 3); % Normalize
end
disp('Projected 2D points (x_projected):');
disp(x_projected);

% Computing the transfer geometric error. Side note: textbook's d(x, y) is the
% Euclidean distance between the inhomogeneous points x and y.
transferError = 0;
for i = 1:size(X, 1)
   difference = x(i, 1:2) - x_projected(i, 1:2);
   transferError = transferError + norm(difference)^2;
end
disp("Transfer Error:" );
disp(transferError);
