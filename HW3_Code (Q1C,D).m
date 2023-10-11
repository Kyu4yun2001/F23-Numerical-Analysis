% Direct Calculation Method
dim = 100;
norm1_t = zeros(size(dim));
norminf_t = zeros(size(dim));

for mult = 1:7
    disp(dim);
    A = randi([0,10], dim, dim);

    tic;
    norm_1 = max(sum(abs(A)));
    norm1_time_elapse = toc;
    norm1_t(mult) = norm1_time_elapse;

    tic;
    norm_inf = max(sum(abs(A),2));
    norminf_time_elapse = toc;
    norminf_t(mult) = norminf_time_elapse;

    dim = dim * 2;
end

% Using MatLab Built In Function
dim_m = 100;
m_norm1_t = zeros(size(dim_m));
m_norminf_t = zeros(size(dim_m));

for mult_m = 1:7
    disp(dim_m);
    A_m = randi([0,10], dim_m, dim_m);

    tic;
    m_norm_1 = norm(A_m,1);
    m_norm1_time_elapse = toc;
    m_norm1_t(mult_m) = m_norm1_time_elapse;

    tic;
    m_norm_inf = norm(A_m,inf);
    m_norminf_time_elapse = toc;
    m_norminf_t(mult_m) = m_norminf_time_elapse;
    
    dim_m = dim_m * 2;
end
    
% Graphing Results
figure;
plot(1:7, norm1_t, 'o-',DisplayName='Norm-1, Manual');
hold on;
plot(1:7, norminf_t,'x-',DisplayName='Norm-inf, Manual');
hold on;
plot(1:7, m_norm1_t, '+-',DisplayName='Norm-1, MatLab');
hold on;
plot(1:7, m_norminf_t,'*-',DisplayName='Norm-inf, MatLab');
title('Norm Values vs. Matrix Size');
xlabel('Matrix size increase iterations');
legend('Location','north');
