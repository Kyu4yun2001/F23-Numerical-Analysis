% 1A
x_1A = -3:0.25:3;
y11_1A = sqrt(1-0.25*x_1A.^2);
y12_1A = -sqrt(1-0.25*x_1A.^2);
y2_1A = (sqrt(3)/2)*x_1A.^2;

plot(x_1A,y11_1A,'b',x_1A,y12_1A, 'b');
hold on;
plot(x_1A,y2_1A,'DisplayName','f2(x,y)');
hold on;
ylim([-2,2]);

% ====================================================================== %

%1D
f = @(x,y) [2*x.^2+8*y.^2-8;y-(sqrt(3)/2)*x.^2];
J_f = @(x,y) [4*x, 16*y;-sqrt(3)*x,1];
num_iter = 5;

% For (x0, y0) = (2,3)
fprintf('For starting points (x0,y0) = (2,3): \n');
x1 = 2;
y1 = 3;

x1_val = zeros(1,num_iter+1);
y1_val = zeros(1,num_iter+1);

for i = 1:5
   x1_val(i) = x1;
   y1_val(i) = y1;

   f_val = f(x1,y1);
   J_val = J_f(x1, y1);

   newton_f_val = [x1;y1] - inv(J_val) * f_val;
   
   x1 = newton_f_val(1,1);
   y1 = newton_f_val(2,1);
   fprintf('Iteration %d: x1 = %f, y1 = %f\n', i, x1, y1);
end

% For (x0, y0) = (-1.5,2)
fprintf('For starting points (x0,y0) = (-1.5,2): \n');
x2 = -1.5;
y2 = 2;

x2_val = zeros(1,num_iter+1);
y2_val = zeros(1,num_iter+1);

for i = 1:5
   x2_val(i) = x2;
   y2_val(i) = y2;

   f_val = f(x2,y2);
   J_val = J_f(x2, y2);

   newton_f_val = [x2;y2] - inv(J_val) * f_val;
   
   x2 = newton_f_val(1,1);
   y2 = newton_f_val(2,1);
   fprintf('Iteration %d: x1 = %f, y1 = %f\n', i, x2, y2);
end

figure;
% from 1A
plot(x_1A,y11_1A,'b',x_1A,y12_1A, 'b');
hold on;
plot(x_1A,y2_1A,'g');
hold on;
% from 1D
plot(x1_val, y1_val, 'm');
hold on;
plot(x2_val,y2_val,'k');

% ====================================================================== %

% Question 5
U1 = [1,2,6,-1;0,3,1,0;0,0,4,-1;0,0,0,2];
b1 = [-1;-3;-2;4];

x = backwards(U1,b1);
disp(x);

function x = backwards(U,b)
    % checking size of Upper Triangle Matrix
    [m_m,m_n] = size(U);
    if m_m ~= m_n
        error("Upper Triangle is not square, please try again");
    end
    
    % Checking size of vector
    [v_n] = size(b);
    if v_n ~= m_n
        error("Matrix U and vector b are not compatible size, please try again");
    end
    
    % Performing backward substitution
    x = U\b;
end


