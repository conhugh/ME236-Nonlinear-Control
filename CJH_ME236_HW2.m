%Connor Hughes
%ME 236 HW 2

close all

%generate a random 10 element vector for theta:
th = 100.*rand(10, 1) - 50;
%generate a random 10 element vector of initial conditions for theta_hat:
th_hat = 100.*rand(10, 1) - 50;
%initialize p to zero:
p = zeros(10, 1);
%initialize parameters:
eps = 0.1;
gam = 100;

%use ODE45 to solve the system for t = 0 -> 200 s:
tspan = [0 200];
t = zeros(20, 1);  %only needed to keep ode45 happy, not used
x0 = [th_hat; p];
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[T, Y] = ode15s(@(t, x) diffeq(t, th, x, eps, gam), tspan, x0, options);
Y = Y';

%calculate V(x, t) = 1/2*(th_hat - th)'*(th_hat - th) + 1/2*(p'*p):
V = zeros(length(Y), 1);
for i = 1:length(Y)
   V(i) = 0.5.*(Y(1:10, i) - th)'*(Y(1:10, i) - th) + 0.5.*Y(11:20, i)'*Y(11:20, i);
end

%plot log10(V(x, t)) vs time
figure
plot(T, log10(V));
title("Gamma = 100 case:")
xlabel("t (seconds)")
ylabel("log10(V)")
hold on

%repeat the process for the gamma = 0 case:
[T, Y] = ode15s(@(t, x) diffeq(t, th, x, eps, 0), tspan, x0, options);
Y = Y';

%calculate V(x, t) = 1/2*(th_hat - th)'*(th_hat - th) + 1/2*(p'*p):
V = zeros(length(Y), 1);
for i = 1:length(Y)
   V(i) = 0.5.*(Y(1:10, i) - th)'*(Y(1:10, i) - th) + 0.5.*Y(11:20, i)'*Y(11:20, i);
end

%plot log10(V(x, t)) vs time
figure
plot(T, log10(V));
title("Gamma = 0 case:")
xlabel("t (seconds)")
ylabel("log10(V)")
hold on

%repeat the process for the system given in part 3b:
th_hat0 = zeros(10, 1);
[T, Y] = ode15s(@(t, th_hat) diffeq2(t, th_hat, th), tspan, th_hat0, options);
Y = Y';

%calculate V(th_hat, th, t) = (th_hat - th)'*(th_hat - th):
V = zeros(length(Y), 1);
for i = 1:length(Y)
   V(i) = (Y(:, i) - th)'*(Y(:, i) - th);
end

%plot log10(V(th_hat, th, t)) vs time
figure
plot(T, log10(V));
title("Part 3b:")
xlabel("t (seconds)")
ylabel("log10(V)")
hold on  


function [x_dot] = diffeq(t, th, x, eps, gam)
    x_dot = zeros(20, 1);
    W = W_(t);
    x_dot(1:10) = W'*W*x(11:20);
    x_dot(11:20) = -eps.*x(11:20) - W'*(W*x(1:10) - W*th) - gam.*(sign((W*x(1:10) - W*th)'*W*x(11:20)) + 1)*x(11:20);
end

function [th_hat_dot] = diffeq2(t, th_hat, th)
    th_hat_dot = zeros(10, 1);
    W = W_(t);
    th_hat_dot = -W'*(W*th_hat - W*th);
end

function [W] = W_(t)
    W = zeros(1, 10);
    W(cast(floor(mod(t, 10)), 'int8') + 1) = 1;
end