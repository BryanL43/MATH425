% Exercise 5
% Part a
A = imread("sillyDog.png");
A = im2gray(A);
A = im2double(A);

% Part b
[P, S, Q] = svd(A);
singularValues = diag(S);
disp('Image Singular values:');
disp(singularValues);

k = 15;
A_k = P(:, 1:k) * S(1:k, 1:k) * Q(:, 1:k)';
imshow(A_k);
