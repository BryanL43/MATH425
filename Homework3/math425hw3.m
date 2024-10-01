% Exercise 1
% 1a
fprintf("Exercise 1a:\n");
v = rand(5, 1);
w = rand(3, 1);
A = v * w';
disp(A);
fprintf("The rank of A is %d\n", rank(A));


% Exercise 2
fprintf("\n");
fprintf("Exercise 2:\n");
v1 = [1; 2; 0; 1];
v2 = [0; -1; 3; 0];
v3 = [2; 0; 1; -1];
v4 = [3; 0; -1; -2];
V = [v1 v2 v3 v4];
disp(V);

[n, k] = size(V);
if (rank(V) == k)
    fprintf("The vectors are linearly independent and not a linear combination of each other.\n");
elseif (rank(V) < k)
   fprintf("The vectors are linear dependent/combination of each other.\n");
end


% Exercise 3
fprintf("\n");
v1 = [1; 0; 2];
v2 = [3; -1; 1];
v3 = [2; -1; -1];
v4 = [4; -1; 3];

fprintf("Exercise 3a:\n");
A = [v1 v2 v3 v4];

if (rank(A) == 3)
    fprintf("The vectors v1, v2, v4, v4 span R^3.\n");
else
    fprintf("The vectors v1, v2, v3, v4 do NOT span R^3.\n");
end
fprintf("Since the rank(A) not equals to 3, we have only a span of R^2 due to v3/v4 being\n" + ...
    "a combination of v1 and v2.\n")

fprintf("\nExercise 3b:\n");
A = [v1 v2 v3 v4];
augmented_A = [A zeros(size(A, 1), 1)];
disp(rref(augmented_A));
fprintf("Here we can see that rref(A) has a free variable, thus not linear independent.\n");
fprintf("Additionally, according to the textbook Lemma 2.23 (pg. 96), it states that if\n" + ...
    "any collection of k (column vectors) > n (row) vectors in R^n then the vectors are linearly dependent.\n");


% Exercise 4
fprintf("\n");
fprintf("Exercise 4:\n");

function B = myGS(A)
    [m, n] = size(A);
    B = zeros(m, n);
    V = A; % Use for intermediate vectors for computing

    % Iterate through matrix columns and perform modified Gram-Schmidt
    for i = 1:n
        % Normalize the current column vector
        B(:, i) = V(:, i) / norm(V(:, i));

        % For each j > i vectors, subtract the projection,
        % i.e. <w_k, u_{k-1}>u_{k-1}
        for j = i+1:n
            V(:, j) = V(:, j) - (B(:, i)' * V(:, j)) * B(:, i);
        end
    end
end

A = [1 0 1 1; 0 1 0 1; 1 0 0 1; 0 -1 1 1];
B = myGS(A);
disp(B);
