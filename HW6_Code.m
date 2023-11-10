% Q3 A & B
Q3_A = [-2 1 4;
    1 1 1;
    4 1 -2];
Iterations = 5;

fprintf("Actual E-Val, E-Vector \n");
[e_vec, e_val] = eig(Q3_A);
disp(e_vec);
disp(e_val);

fprintf('For x0 = [1, 2, -1]^T \n');
x0 = [1; 2; -1];
[x_k1_iter, y_iter] = powerMethod(Q3_A, x0, Iterations);

fprintf('For x0 = [1, 2, 1]^T \n');
x0 = [1; 2; 1];
[x_k1_iter, y_iter] = powerMethod(Q3_A, x0, Iterations);

% Q3 C & D
fprintf("For Inital v_guess = -5\n")
v_guess = -5;
x0 = [1; 0;0];
[y_k] = inverseIteration(Q3_A, x0, v_guess, Iterations);
disp(y_k);

fprintf("For Initial v_guess = 1\n")
v_guess = 1;
x0 = [1; 0;0];
[y_k] = inverseIteration(Q3_A, x0, v_guess, Iterations);
disp(y_k);

fprintf("For Initial v_guess = 4 \n")
v_guess = 4;
x0 = [1; 0;0];
[y_k] = inverseIteration(Q3_A, x0, v_guess, Iterations);
disp(y_k);

% Q4 - Tridiagonalization w Householders
Q4_A = [2 1 2 2;
        1 -7 6 5;
        2 6 2 -5;
        2 5 -5 1];

[T, Q] = tridiagonalize(Q4_A);

fprintf('Tridiagonal Matrix T:\n');
disp(T);

% Functions

% Q3A - Power Method Function
function [x_k1_iter, y_iter] = powerMethod(A, x0, Iterations)
    x = x0;
    x_k1_iter = zeros(length(x0), Iterations);
    y_iter = zeros(length(x0), Iterations);

    for k = 1:Iterations
        x_k1 = A * x;
        y_k = x_k1 / norm(x_k1);

        x_k1_iter(:,k) = x_k1;
        y_iter(:,k) = y_k;

        x = x_k1;
    end

    disp(x_k1_iter);
    disp(y_iter);
end

% Q3C - Inverse Iteration Method
function [inverse_y_iter] = inverseIteration(A, x0, v_guess, Iterations)
    x = x0;
    inverse_y_iter = zeros(length(x0), Iterations);

    for k = 1:Iterations
        n = size(A, 1);
        I = eye(n);

        [L, U, P] = lu(A - v_guess * I);
        x_k1 = U \ (L \ (P * x));
        y_k = x_k1 / norm(x_k1);

        inverse_y_iter(:,k) = y_k;

        x = y_k;
    end 
end


% Q4 - Tridiagonalization Function
function [T, Q] = tridiagonalize(A)
    [m, n] = size(A);
    T = A;
    Q = eye(m);
    
    for k = 1:(n-2)
        x = T((k+1):end, k);
        v = zeros(size(x));
        v(1) = sign(x(1)) * norm(x) + x(1);
        v(2:end) = x(2:end);
        v = v / norm(v);
        
        H = eye(m);
        H((k+1):end, (k+1):end) = H((k+1):end, (k+1):end) - 2 * v * v';
        
        T = H * T * H;
        Q = Q * H';
    end
end

