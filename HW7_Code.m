%% 
% 1B - Simple QR iteration
r = 0.25;
A = [1, r;r, 1];
Ak = A;
numOfIterations = 1e5; 
for i=1:numOfIterations
    [Q,R] = qr(Ak);
    Ak = R*Q;
end

disp(Ak)

%%
% 1C - QR Algorithm with Tolerance
tolerance = 1e-10;
[eigenvalues, eigenvectors, iterations] = qr_tol(A, tolerance   );

disp('Eigenvalues:');
disp(eigenvalues);
disp('Eigenvectors:');
disp(eigenvectors);
disp(['Number of iterations: ', num2str(iterations)]);

%%
% 1D - QR Algorithm for multiple values of r
r = [0.1, 0.3, 0.5, 0.7, 0.9];
iter_vec = zeros(5);

for i = 1:length(r)
    r_val = r(i);
    disp(r_val);
    Q1D_A = [1, r_val;r_val, 1];
    [e_val, e_vec, iter] = qr_tol(Q1D_A, 1e-10);
    iter_vec(i) = iter;
end    

figure;
plot(r, iter_vec, '-o');
xlabel('r');
ylabel('Number of Iterations');
grid on;

%%
function [qr_evalues, qr_evectors, qr_iter] = qr_tol(A, tolerance)
    n = size(A, 1);
    Q = eye(n);
    Ak = A;
    max_iterations = 1e5;

    for qr_iter = 1:max_iterations
        [Q, R] = qr(Ak);
        Ak = R*Q;
        if max(abs(triu(Ak, 1))) < tolerance
            break;
        end
    end

    if qr_iter == max_iterations
        error('QR algorithm did not converge.');
    end

    qr_evalues = diag(Ak);
    qr_evectors = Q;
end
