%Connor Hughes
%ME 236 HW 4

%Problem 10.14:

%initialize TORA parameters:
M = 1.3608;
m = 0.096;
L = 0.0592;
J = 0.0002175;
k = 186.3;

%set initial conditions for TORA state
x0 = [pi; 0; 0.025; 0];

%% Passivity-Based Control
%initialize control parameters
U1 = 0.07;
U2 = 0.1 - U1;
k1 = 1.0;
k2 = 0.085;

%use ODE45 to solve the system for t = 0 -> 200 s with Passivity-Based Control:
tspan = [0 15];
t = zeros(4, 1);  %only needed to keep ode45 happy, not used
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[T, Y] = ode15s(@(t, x) TORAPassive(x, m, M, L, J, k, k1, k2, U1, U2), tspan, x0, options);
Y = Y';

% %plot components of x vs time
% figure
% plot(T, Y(1, :));
% title("Passivity-Based Control:")
% xlabel("t (seconds)")
% ylabel("$\theta$", 'Interpreter', 'latex')
% hold on
% 
% figure
% plot(T, Y(2, :));
% title("Passivity-Based Control:")
% xlabel("t (seconds)")
% ylabel("$\dot{\theta}$", 'Interpreter', 'latex')
% 
% figure
% plot(T, Y(3, :));
% title("Passivity-Based Control:")
% xlabel("t (seconds)")
% ylabel("$x_c$", 'Interpreter', 'latex')
% 
% figure
% plot(T, Y(4, :));
% title("Passivity-Based Control:")
% xlabel("t (seconds)")
% ylabel("$\dot{x_c}$", 'Interpreter', 'latex')

%% Sliding Mode Control
%initialize control parameters:
beta = 0.1;
k1 = 5;
k2 = 1100;
mu = 10;

%use ODE15s to solve the system for t = 0 -> 200 s with Sliding Mode Control:
tspan = [0 15];
t = zeros(4, 1);  %only needed to keep ode45 happy, not used
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[T, Y] = ode15s(@(t, x) TORASMC(x, m, M, L, J, k, beta, k1, k2, mu), tspan, x0, options);
Y = Y';

% %plot components of x vs time
% figure
% plot(T, Y(1, :));
% title("Sliding Mode Control:")
% xlabel("t (seconds)")
% ylabel("$\theta$", 'Interpreter', 'latex')
% 
% figure
% plot(T, Y(2, :));
% title("Sliding Mode Control:")
% xlabel("t (seconds)")
% ylabel("$\dot{\theta}$", 'Interpreter', 'latex')
% 
% figure
% plot(T, Y(3, :));
% title("Sliding Mode Control:")
% xlabel("t (seconds)")
% ylabel("$x_c$", 'Interpreter', 'latex')
% 
% figure
% plot(T, Y(4, :));
% title("Sliding Mode Control:")
% xlabel("t (seconds)")
% ylabel("$\dot{x_c}$", 'Interpreter', 'latex')

function [u] = uPassivity(x1, x2, k1, k2, U1, U2) 
    u = -U1*min(max(-1, k1*x1), 1) - U2*min(max(-1, k2*x2), 1);
end

function [u] = uSMC(beta, s, mu) 
    u = -beta*min(max(-1, s/mu), 1);
end

function [del] = delta(x1, m, M, L, J)
    del = (J + m*L^2)*(m + M) - m^2*L^2*(cos(x1))^2;
end

function [x_dot] = TORAPassive(x, m, M, L, J, k, k1, k2, U1, U2)
    x_dot = zeros(4, 1);
    x_dot(1) = x(2);
    x_dot(2) = (1/delta(x(1), m, M, L, J))*((m + M)*(uPassivity(x(1), x(2), k1, k2, U1, U2)) - m*L*cos(x(1))*(m*L*x(2)^2*sin(x(1))-k*x(3)));
    x_dot(3) = x(4);
    x_dot(4) = (1/delta(x(1), m, M, L, J))*(-m*L*(uPassivity(x(1), x(2), k1, k2, U1, U2))*cos(x(1)) + (J + m*L^2)*(m*L*x(2)^2*sin(x(1)) - k*x(3)));
end

function [x_dot] = TORASMC(x, m, M, L, J, k, beta, k1, k2, mu)
    s = x(2) + k1*x(1) - k2*x(3)*cos(x(1));
    x_dot = zeros(4, 1);
    x_dot(1) = x(2);
    x_dot(2) = (1/delta(x(1), m, M, L, J))*((m + M)*(uSMC(beta, s, mu)) - m*L*cos(x(1))*(m*L*x(2)^2*sin(x(1))-k*x(3)));
    x_dot(3) = x(4);
    x_dot(4) = (1/delta(x(1), m, M, L, J))*(-m*L*(uSMC(beta, s, mu))*cos(x(1)) + (J + m*L^2)*(m*L*x(2)^2*sin(x(1)) - k*x(3)));
end
