%The Newton Raphson Method
clear, clc
%input Sections
a = 0.1 ; %OA=0.1m
b = 0.2 ; %AB=0.2m
o=1; %omega=1rad/s
N = 101 ; %number of steps 
u0 = [pi ; 0]; %starting point
i = 1;
u = u0;
err = 1e-9;

%create matrix
theta = zeros(N, 1); %angle 
x_matrix = zeros(N, 1); %displacement
theta_derivativ = zeros(N, 1); %angle derivative
x_derivativ = zeros(N, 1); %displacement derivative

for t = linspace(0, 1, N)
    phi = pi / 6 + o * t;
    F = @(u) constraint(u, a, b, phi);
    J = @(u) jacobian(u, b);
    
    [u, counter] = NR_method(F, J, u, err); 
     
    df_dt = [- a * o * sin(o * t)
            a * o * cos(o * t)];
    df_dq = [- b * sin(u(1)), -1
            - b * cos(u(1)), 0];
    u_derivative = - df_dt \ df_dq;
    
    theta(i) = u(1);
    x_matrix(i) = u(2);
    theta_derivativ(i) = u_derivative(1);
    x_derivativ(i) = u_derivative(2);
    
    i = i + 1;
end

%plotting
t = linspace(0, 1, N);
subplot(2,2,1)
plot(t, theta);
title('𝑡 versus angle theta');
subplot(2,2,2)
plot(t, x_matrix);
title('𝑡 versus displacement');
subplot(2,2,3)
plot(t, theta_derivativ);
title('𝑡 versus theta derivative');
subplot(2,2,4)
plot(t, x_derivativ);
title('𝑡 versus displacement derivative');

%functions
function P = constraint(u, a, b, phi)
    theta = u(1);
    d = u(2);
    P = [a * cos(phi) + b * cos(theta) - d
    a * sin(phi) - b * sin(theta)];
end
function P = jacobian(u, b)
    theta = u(1);
    P=[- b * sin(theta), -1
        - b * cos(theta), 0];
end
%Newton-Raphson’s method
function [x, counter] = NR_method(F, J, x, eps)
    F_value = F(x);
    F_norm = norm(F_value); % L2 norm of vector
    counter = 0;
    while F_norm > eps && counter < 100
        delta = J(x) \ - F_value;
        x = x + delta;
        F_value = F(x);
        F_norm = norm(F_value);
        counter = counter + 1;
    end
    if F_norm > eps
        counter = -1;
    end
end