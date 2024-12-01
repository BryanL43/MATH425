% Exercise #6
function [Q_, E_] = computeEigenValues(A)
    E_ = A; % Make a copy of A

    Q_hat = eye(size(E_));
    
    % Repeatedly perform QR-factorization to compute the eigenvalues
    while true
        [Q, R] = qr(E_);
        Q_hat = Q_hat * Q;
        E_ = R * Q;
        
        % End loop when non-diagonal entries are less than the defined
        % error margin
        offDiagonalEntries = E_ - diag(diag(E_));
        if all(abs(offDiagonalEntries(:)) < 1e-6)
            break;
        end
    end

    Q_ = Q_hat;
end

% Test #1
disp("-------------------- [Test #1] --------------------");
A = [3 1 2; 1 3 1; 2 1 3];
[Q, L] = eig(A);
disp("MATLAB computed eigenvectors");
disp(Q);
disp("MATLAB computed eigenvalues:");
disp(L);

[Q_, L_] = computeEigenValues(A);
disp("QR Algorithm eigenvectors:");
disp(Q_);
disp("QR Algorithm eigenvalues:");
disp(L_);

% Test #2
disp("-------------------- [Test #2] --------------------");
A = [31 -1 30 -9; -1 14 -2 -1; 30 -2 31 -4; -9 -1 -4 22];
[Q, L] = eig(A);
disp("MATLAB computed eigenvectors");
disp(Q);
disp("MATLAB computed eigenvalues:");
disp(L);

[Q_, L_] = computeEigenValues(A);
disp("QR Algorithm eigenvectors:");
disp(Q_);
disp("QR Algorithm eigenvalues:");
disp(L_);
