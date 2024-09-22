% Exercise 1
function U = myPartialPivot(A)
    rowInterchanged = 0;
    [n, m] = size(A);
    if n ~= m
        error("Expected A to be an nxn matrix.");
    end

    U = A;

    for i = 1:n-1
        [~, maxIndex] = max(abs(U(i:n, i))); % row index of largest absolute element in column

        % Get row index of largest absolute element relative to U.
        % The previous statement is a submatrix of U, thus need to require
        % relative.
        maxIndex = maxIndex + i - 1;

        % Swap rows if the pivot is not the largest
        if maxIndex ~= i
            U([i, maxIndex], :) = U([maxIndex, i], :);
            rowInterchanged = rowInterchanged + 1;
        end
        
        % Perform Gaussian elimination below the pivot
        for j = i + 1:n
            factor = U(j, i) / U(i, i);
            U(j, i:n) = U(j, i:n) - factor * U(i, i:n);
        end
    end
    fprintf("Number of row interchanges when doing partial pivoting: %d\n", rowInterchanged);
end

function rank = myRank(A)
    [n, m] = size(A);
    if n ~= m
        error("Expected A to be an nxn matrix.");
    end

    U = myPartialPivot(A);
    rank = 0;
    
    % Iterate through each row of the upper-triangular matrix and counting
    % non-zero rows to acquire rank of the matrix.
    for i = 1:n
        if any(round(U(i, :), 10) ~= 0) % Using round to avoid floating-point precision issues
            rank = rank + 1;
        end
    end
end

P = rand(5, 3);
Q = rand(3, 5);
A = P * Q;
fprintf("Rank of A = P * Q: %d\n\n", myRank(A));


% Exercise 2
disp("Example of strictly diagonal dominant matrix, which is not a diagonal matrix:");
A = [5 2 1 0; 1 6 1 2; 0 2 7 3; 0 0 1 8];
disp(A);
disp("We can see that there are non-zero elements outside the diagonals, thus not diagonal matrix.");
fprintf("\n");

fprintf("No row interchanges with diagonally dominant matrix (shown below by partial pivoting the matrix A above):\n")
myPartialPivot(A);

fprintf("\nAnother example to satisfy 2b:\n");
C = [5 1 1; 1 6 1; 1 1 7];
disp(C);
myPartialPivot(C);
