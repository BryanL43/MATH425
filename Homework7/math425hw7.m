% Exercise 2 (further justification via example)
% A is a full-rank matrix
A = [1 2 2;
     2 1 2;
     2 2 1];

[~, S, ~] = svd(A);
disp("MATLAB svd acquired singular values:");
disp("Regular S singular values:");
disp(S);
disp("Inverse S singular values:");
disp(inv(S));

K = A' * A;
[Q1_eig, L1_eig] = eig(K);

K_1 = inv(A)' * inv(A);
[Q2_eig, L2_eig] = eig(K_1);

disp("Computed K = A' * A eigenvalue:");
disp(L1_eig);
disp("Computed K associated singular value:");
disp(sqrt(L1_eig));
disp("Computed K_1 = inv(A)' * inv(A) eigenvalue:");
disp(L2_eig);
disp("Computed K_1 associated singular value:");
disp(sqrt(L2_eig));

% Exercise 5
A = imread("sillyDog.png");
A = im2gray(A);
A = im2double(A);

[P, S, Q] = svd(A);
singularValues = diag(S);
disp('Image Singular values:');
disp(singularValues);

k = 15;
A_k = P(:, 1:k) * S(1:k, 1:k) * Q(:, 1:k)';
imshow(A_k);
