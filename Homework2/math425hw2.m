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
    disp("Number of row interchanges when doing partial pivoting:")
    disp(rowInterchanged);
end

function rank = myRank(A)
    rank = 0;
    U = myPartialPivot(A);
    [n, m] = size(U);
    
    for i = 1:n
        if norm(U(i, :)) > 1e-12
            rank = rank + 1;
        end
    end
end

%A = [1 0.6; 0.01 1.6];
%A = [0.02 0.01 0 0; 1 2 1 0; 0 1 2 1; 0 0 100 200];

P = rand(5, 3);
Q = rand(3, 5);
A = P * Q;
disp(myRank(A));


% Exercise 2
disp("Example of strictly diagonal dominant matrix, which is not a diagonal matrix:");
A = [5 2 1 0; 1 6 1 2; 0 2 7 3; 0 0 1 8];
disp(A);

%C = [-4 2 1; 1 6 2; 1 -2 5];
%disp(C);
%myPartialPivot(C);

myPartialPivot(A);
