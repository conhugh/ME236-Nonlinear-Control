%Connor Hughes
%ME 236 HW 3

close all

%initialize TORA parameters
M = 1.3608;
m = 0.096;
L = 1;
J = 0.0002175;
k = 186.3;

%set initial conditions for TORA state
x0 = rand(4, 1).*10 - 5;
x02 = rand(5, 1).*10 - 5;

%use ODE45 to solve the system for t = 0 -> 200 s:
tspan = [0 200];
t = zeros(4, 1);  %only needed to keep ode45 happy, not used
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[T, Y] = ode45(@(t, x) TORA(x, -x(2), m, M, L, J, k), tspan, x0, options);
Y = Y';

%calculate V(x)
V = zeros(length(Y), 1);
for i = 1:length(Y)
   V(i) = Lyap(Y(:, i), m, M, L, J, k);
end

%plot log10(V(x)) vs time
figure
%%
% 
%   for x = 1:10
%       disp(x)
%   end
% 
plot(T, log10(V));
title("Feedback: w = -x2:")
xlabel("t (seconds)")
ylabel("log10(V)")
hold on

% %plot xc vs theta
figure
plot(Y(1, :), Y(2, :));
title("Feedback: w = -x2:")
xlabel("theta")
ylabel("xc")
hold on

%now for the feedback dynamics of problem 4b
%use ODE45 to solve the system for t = 0 -> 500 s:
tspan = [0 500];
t = zeros(5, 1);  %only needed to keep ode45 happy, not used
[T, Y] = ode45(@(t, x) TORA2(x, m, M, L, J, k), tspan, x02, options);
Y = Y';

%calculate V(x)
V = zeros(length(Y), 1);
for i = 1:length(Y)
   V(i) = Lyap2(Y(:, i), m, M, L, J, k);
end

%plot log10(V(x)) vs time
figure
plot(T, log10(V));
title("Problem 4 feedback dynamics:")
xlabel("t (seconds)")
ylabel("log10(V)")
hold on

%plot xc vs theta
figure
plot(Y(1, :), Y(2, :));
title("Problem 4 feedback dynamics:")
xlabel("theta")
ylabel("xc")
hold on


%will use phi_1(x1) = x1
function [x_dot] = TORA(x, w, m, M, L, J, k)
    x_dot = zeros(4, 1);
    x_dot(1) = x(2);
    x_dot(2) = (1/delta(x(1), m, M, L, J))*((m + M)*(-x(1) + w) - m*L*cos(x(1))*(m*L*x(2)^2*sin(x(1))-k*x(3)));
    x_dot(3) = x(4);
    x_dot(4) = (1/delta(x(1), m, M, L, J))*(-m*L*(-x(1) + w)*cos(x(1)) + (J + m*L^2)*(m*L*x(2)^2*sin(x(1)) - k*x(3)));
end

%will use phi_2(x2)= x2
function [x_dot] = TORA2(x, m, M, L, J, k)
    x_dot = zeros(4, 1);
    x_dot(1) = x(2);
    x_dot(2) = (1/delta(x(1), m, M, L, J))*((m + M)*(-x(1) + x(5)) - m*L*cos(x(1))*(m*L*x(2)^2*sin(x(1))-k*x(3)));
    x_dot(3) = x(4);
    x_dot(4) = (1/delta(x(1), m, M, L, J))*(-m*L*(-x(1) + x(5))*cos(x(1)) + (J + m*L^2)*(m*L*x(2)^2*sin(x(1)) - k*x(3)));
    x_dot(5) = -x(5) - x(2);
end

function [del] = delta(x1, m, M, L, J)
    del = (J + m*L^2)*(m + M) - m^2*L^2*(cos(x1))^2;
end

function [V] = Lyap(x, m, M, L, J, k)
    V = 0.5*(J + m*L^2)*x(2)^2 + m*L*cos(x(1))*x(2)*x(4) + 0.5*(m + M)*x(4)^2 + 0.5*k*x(3)^2 + 0.5*x(1)^2;
end

function [V] = Lyap2(x, m, M, L, J, k)
    V = 0.5*(J + m*L^2)*x(2)^2 + m*L*cos(x(1))*x(2)*x(4) + 0.5*(m + M)*x(4)^2 + 0.5*k*x(3)^2 + 0.5*x(1)^2 + 0.5*x(2)^2;
end