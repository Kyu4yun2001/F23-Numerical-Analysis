% Q1

A1 = [9, -6; 12, -8; 0, 20];
[Q1, R1] = qr(A1);
disp(Q1);
disp(R1);
Q1_hat = Q1(1:3,1:2);
disp(Q1_hat);
R1_hat = R1(1:2,1:2);
disp(R1_hat);

b1 = [300; 600; 900];
x1 = R1_hat \ (Q1_hat' * b1);

x1_var = 50:100;
y11 = (300-9*x1_var)/(-6);
y12 = (600 - 12*x1_var)/(-8);
y13 = (900/20)*ones(size(x1_var));

figure;
plot(x1_var, y11);
hold on;
plot(x1_var, y12);
hold on;
plot(x1_var, y13);
hold on;
plot(x1(1,1), x1(2,1),'ro','MarkerSize',10);

% Q2
A2 = [exp(0), 0.0^2, 0.0, 1; 
    exp(0.5), 0.5^2, 0.5, 1;
    exp(1.0), 1.0^2, 1.0, 1;
    exp(1.5), 1.5^2, 1.5, 1;
    exp(2.0), 2.0^2, 2.0, 1;
    exp(2.0), 2.0^2, 2.0, 1;
    exp(2.5), 2.5^2, 2.5, 1;
    ];
b2 = [0; 0.2; 0.27; 0.3; 0.32; 0.35; 0.27];

x2 = inv(A2' * A2) * A2' * b2;

x2_var = -0.5:4;
y2 = x2(1,1)*exp(x2_var) + x2(2,1)*(x2_var.^2) + x2(3,1)*x2_var + x2(4,1);

figure;
plot(x2_var, y2);
hold on;
plot(0.0, 0.0, 'ro');
hold on;
plot(0.5, 0.2, 'ro');
hold on;
plot(1.0, 0.27, 'ro');
hold on;
plot(1.5, 0.3, 'ro');
hold on;
plot(2.0, 0.32, 'ro');
hold on;
plot(2.0, 0.35, 'ro');
hold on;
plot(2.5, 0.27, 'ro');
hold on;

