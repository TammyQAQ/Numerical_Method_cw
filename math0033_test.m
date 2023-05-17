%% Yitian (Tammy) Wang
% MATH0033 Numerical Methods Computational homework 1
clear all, close all, clc, format long, format compact
%%
% *Exercise 1*
%
% (a)

% Define the function and exact solution
f = @(t, y) 1 - y^2;
yExact = @(t) tanh(t);

% Set the time interval and the number of time steps
T = 20;
N = [19, 21, 40, 80, 160];

% Set the initial condition
y0 = 0;

% Compute the numerical solutions using the forward Euler method
for i = 1:length(N)
    % Set the time step size h
    h = T/N(i);
    
    % Initialize the arrays to store the time points and the numerical solution
    t = zeros(1, N(i)+1);
    y = zeros(1, N(i)+1);
    
    % Set the initial time and solution
    t(1) = 0;
    y(1) = y0;
    
    % Iterate the forward Euler method
    for j = 1:N(i)
        t(j+1) = t(j) + h;
        y(j+1) = y(j) + h*f(t(j), y(j));
    end
    
    % Plot the numerical solution and the exact solution
    figure;
    plot(t, y, t, yExact(t));
    legend('Numerical', 'Exact');
    title(sprintf('N = %d', N(i)));
end

%%
% Compute the errors
efe = zeros(1, N(i)+1);
for j = 1:N(i)+1
    efe(j) = abs(y(j) - yExact(t(j)));
end
e_fe_end(i) = abs(y(end) - yExact(T));

%%
% (b)

% Compute the numerical solutions using Heun's method
for i = 1:length(N)
    % Set the time step size h
    h = T/N(i);
    
    t = zeros(1, N(i)+1);
    y = zeros(1, N(i)+1);
    t(1) = 0;
    y(1) = y0;
    
    % Iterate Heun's method
    for j = 1:N(i)
        t(j+1) = t(j) + h;
        k1 = h*f(t(j), y(j));
        k2 = h*f(t(j) + h, y(j) + k1);
        y(j+1) = y(j) + (k1 + k2)/2;
    end
    
    % Compute the errors
    eheun = zeros(1, N(i)+1);
    for j = 1:N(i)+1
        eheun(j) = abs(y(j) - yExact(t(j)));
    end
    eheun_end(i) = abs(y(end) - yExact(T));
    
    % Plot the numerical solution and the exact solution
    figure;
    plot(t, y, t, yExact(t));
    legend('Numerical', 'Exact');
    title(sprintf('N = %d', N(i)));
end

