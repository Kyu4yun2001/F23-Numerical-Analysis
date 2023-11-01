% HW 5 Question 3A
matrix_size = [5, 10, 20];

for n = matrix_size
    fprintf('Matrix Size:%.1f\n',n);

    H = hilb(n);
    e_n = ones(n, 1);
    norm_1 = norm(H*(inv(H)*e_n) - e_n,1);
    norm_2 = norm(H*(inv(H)*e_n) - e_n,2);
    norm_inf = norm(H*(inv(H)*e_n) - e_n,inf);
    
    fprintf('1-norm: %.20f\n', norm_1);
    fprintf('2-norm: %.20f\n', norm_2);
    fprintf('inf-norm: %.20f\n', norm_inf);
end
%% 

% HW 5 Question 3B

k_values = zeros(1, 16);
derivative_approximate_values = zeros(1, 16);
finite_diff = zeros(1,16);

for k = 1:16
    delta_x = 10^(-k);
    approximate_deri = (my_function(pi/4 + delta_x) - my_function(pi/4)) / delta_x;
    
    k_values(k) = delta_x;
    derivative_approximate_values(k) = approximate_deri;
    finite_diff(k) = abs(approximate_deri - 3.101766393836051);
end

loglog(k_values,finite_diff,'-o');
xlabel('log(k)');
ylabel('log(finite diff = approximate - actual)')

function result = my_function(x)
    numerator = exp(x);
    denominator = (cos(x))^3 + (sin(x))^3;
    result = numerator / denominator;
end




