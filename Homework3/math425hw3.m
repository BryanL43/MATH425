% Exercise 2
fprintf("\n");
fprintf("Exercise 2:\n");
v1 = [1; 2; 0; 1];
v2 = [0; -1; 3; 0];
v3 = [2; 0; 1; -1];
b = [3; 0; -1; -2];
A = [v1 v2 v3 b];
disp(rref(A));
disp("Taking the rref of the augmented matrix A = [ v1 v2 v3 | b ], we can see that the last row has 0s beside the b column.");
disp("This means that the matrix is inconsistent and b is not a linear combination of v1, v2, and v3.");


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
fprintf("Since the rank of A does not equals to 3, but instead it is %d.\n" + ...
    "the vectors spans R^%d.\n", rank(A), rank(A));
disp(rref(A));
fprintf("Taking the rref of matrix A also shows that v3 and v4 are a combination of\n" + ...
    "v1 and v2, where v3 = (-1 * v1) + v2 and v4 = v1 + v2.\n");

fprintf("\nExercise 3b:\n");
disp(rref(A));
fprintf("Here we can see that rref(A) has a free variable (shown by the last row with all 0s), thus not linear independent.\n");
fprintf("Additionally, according to the textbook Lemma 2.23 (pg. 96), it states that if\n" + ...
    "any collection of k (column vectors) > n (row) vectors in R^n then the vectors are linearly dependent.\n");

fprintf("\nExercise 3c:\n");
fprintf("No, v1, v2, v3, v4 do not form a basis for R^3. Since the vectors only span R^2 as shown in 3a\n" + ...
    "and aren't linearly independent as shown in 3b, then we cannot form a basis of R^3.\n");
fprintf("Additionally, it is also not possible to choose some subset which is a basis, as the span of the\n" + ...
    "vectors is only R^2.\n");

fprintf("\nExercise 3d:\n");
fprintf("The dimension of the span of v1, v2, v3, v4 is 2. As explained in 3a, there are only 2 linearly independent\n" + ...
    "vectors in the set, shown by the rank(A) = 2.\n");


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
