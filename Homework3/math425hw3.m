% Exercise 2
fprintf("\n");
fprintf("Exercise 2:\n");
v1 = [1; 2; 0; 1];
v2 = [0; -1; 3; 0];
v3 = [2; 0; 1; -1];
b = [3; 0; -1; -2];
A = [v1 v2 v3 b];
disp(rref(A));
fprintf("Taking the rref of the augmented matrix A = [ v1 v2 v3 | b ],\n" + ...
    "we can see that the last row has 0s beside the b column, i.e. 0 = 1.\n" + ...
    "This means that the matrix is inconsistent and b is not a linear combination of v1, v2, and v3.\n");


% Exercise 3
fprintf("\n");
v1 = [1; 0; 2];
v2 = [3; -1; 1];
v3 = [2; -1; -1];
v4 = [4; -1; 3];

fprintf("Exercise 3a:\n");
A = [v1 v2 v3 v4];
fprintf("To span R^3, we need a set of three linearly indendent vectors\n" + ...
    "which can be identified by determining the rank.\n");
disp("rref of A:")
disp(rref(A));
fprintf("The rank of the augmented matrix A has a rank of %d.\n", rank(A));

fprintf("Since the rank of A does not equals to 3, but instead it is %d.\n" + ...
    "The vectors do not spans R^3.\n", rank(A));

fprintf("\nExercise 3b:\n");
disp(rref(A));
fprintf("Here we can see that rref(A) has a free variable (shown by the last row with all 0s), thus not linear independent.\n");
fprintf("Additionally, according to the textbook Lemma 2.23 (pg. 96), it states that if\n" + ...
    "any collection of k (column vectors) > n (row) vectors in R^n then the vectors are linearly dependent.\n");

fprintf("\nExercise 3c:\n");
fprintf("No, v1, v2, v3, v4 do not form a basis for R^3 because the vectors are linearly dependent.\n" + ...
    "Shown in 3a, we found that the matrix comprised of the column vectors has a rank of 2.\n" + ...
    "This means that there are only two linearly independent vectors in the set and are not sufficient to span R^3\n" + ...
    "nor any subset.\n");

fprintf("\nExercise 3d:\n");
fprintf("The dimension of the span of v1, v2, v3, and v4 is 2. Since the dimension reflects the number of vectors in the basis,\n" + ...
    "we can determine this by examining the rank of the augmented matrix A. As demonstrated in part 3a, the rank is 2,\n" + ...
    "which indicates that the basis of the subspace spanned by these vectors must contain 2 vectors.\n" + ...
    "Therefore, we conclude that dim(span{v1, v2, v3, v4}) = 2.\n");


% Exercise 4
fprintf("\n");
fprintf("Exercise 4:\n");

function B = myGS(A)
    [m, n] = size(A);
    B = zeros(m, n);
    V = zeros(m, n);

    % Iterate through augmented A matrix columns to compute orthogonal set
    % {v_1, ..., v_n}
    for i = 1:n
        v = A(:, i); % Acquire i-th column of A
        
        % For 2nd column vector and beyond, perform Gram-Schmidt process
        for j = 1:i-1 % j loop will acquire previous vector
            v = v - ( ( dot(v, V(:, j)) / dot(V(:, j), V(:, j)) ) * V(:, j) );
        end

        V(:, i) = v; % Temp store the v to normalize at the end
    end
    
    % Normalize to acquire orthonormal basis
    for i = 1:n
        B(:, i) = V(:, i) / norm(V(:, i));
    end
end

function B = myGS2(A)
    [m, n] = size(A);
    B = zeros(m, n);

    % Iterate through augmented A matrix columns
    for i = 1:n
        v = A(:, i); % Acquire i-th column of A
        
        % Orthogonalize v against all previous orthonormal vectors
        for j = 1:i-1
            v = v - dot(v, B(:, j)) * B(:, j);
        end
        
        % Normalize to acquire orthonormal vector
        B(:, i) = v / norm(v);
    end
end

w1 = [1; 0; 1; 0];
w2 = [0; 1; 0; -1];
w3 = [1; 0; 0; 1];
w4 = [1; 1; 1; 1];

A = [w1 w2 w3 w4];
B = myGS(A);
disp("Regular GS:")
disp(B);

A2 = [w1 w2 w3 w4];
B2 = myGS2(A2);
disp("Modified GS:")
disp(B2);
