%Connor Hughes
%ME 236 
%Homework 1

close all

%generate a random 20x20 symmetric, positive-definite matrix:
R = rand(20);
P = transpose(R)*R;
%generate a random 20x1 vector of initial conditions for p:
p0 = 100.*rand(20, 1) - 50;
%generate initial conditions for q:
q0 = zeros(20, 1);
%initialize parameters:
eps = 0.1;
mu = 100;

%use ODE45 to solve the system for t = 0 -> 200 s:
tspan = [0 200];
t = zeros(20, 1);  %only needed to keep ode45 happy, not used
x0 = [p0; q0]; 
[T, Y] = ode45(@(t, x) diffeq(t, x, P, eps, mu), tspan, x0);
Y = transpose(Y);
p = Y(1:20, :);
q = Y(21:40, :);

%calculate V(p, q) =(1/2)q_transpose*q + phi(p)
%                  =(1/2)q_transpose*q + p_transpose*P*p:
V = zeros(length(p), 1);
for i = 1:length(p)
   V(i) = 0.5.*transpose(q(:, i))*q(:,i) + transpose(p(:, i))*P*p(:, i);
end

%plot log10(V(p, q)) vs time
figure
plot(T, log10(V));
title("Mu = 100 Case")
xlabel("t (seconds)")
ylabel("log10(V)")
hold on

%set mu = 0, repeat the above:
mu = 0;
tspan = [0 200];
t = zeros(20, 1);  %only needed to keep ode45 happy, not used
[T, Y] = ode45(@(t, x) diffeq(t, x, P, eps, mu), tspan, x0);
Y = transpose(Y);
p = Y(1:20, :);
q = Y(21:40, :);

V = zeros(length(p), 1);
for i = 1:length(p)
   V(i) = 0.5.*transpose(q(:, i))*q(:,i) + transpose(p(:, i))*P*p(:, i);
end

figure
plot(T, log10(V));
title("Mu = 0 Case")
xlabel("t (seconds)")
ylabel("log10(V)")

%simulate and plot the Lyapunov function V = transpose(p)*P*p
%for the system p_dot = -grad_phi(p):
tspan = [0 200];
t = zeros(20, 1); 
[T, Y] = ode45(@(t, x) diffeq2(t, x, P), tspan, p0);
p = transpose(Y);

V = zeros(length(p), 1);
for i = 1:length(p)
   V(i) = transpose(p(:, i))*P*p(:, i);
end

figure
plot(T, log10(V));
xlabel("t (seconds)")
ylabel("log10(V)")


function [x_dot] = diffeq(t, x, P, eps, mu)
    x_dot = zeros(40, 1);
    x_dot(1:20) = x(21:40);
    x_dot(21:40) = -eps.*x(21:40) - 2.*P*x(1:20) - mu.*(sign(transpose(2.*P*x(1:20))*x(21:40)) + 1).*x(21:40);
end

function [x_dot] = diffeq2(t, x, P) 
    x_dot = -2.*P*x(1:20);
end

