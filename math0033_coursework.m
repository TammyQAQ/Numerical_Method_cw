%% Yitian (Tammy) Wang
% MATH0033 Numerical Methods Computational homework 1
clear all, close all, clc, format long, format compact
%%
% *Exercise 1*
%
% (a)
%
% Use Bisection method and plot on the interval between -pi to pi.

x = linspace (-pi, pi);
f = @(x) x/2-sin(x)+pi/6-sqrt(3)/2;
plot(x, f(x))
[zero, res, niter, itersb] = bisection(f, pi/2, pi, 1e-10, 100)

%%
% (b)
%
% From the previous worksheet, the positive root has quadratic convergence
% where the negative root has linear convergence. So the negative root
% would need more number of iterations.


df=@(x) 1/2 - cos(x);
[zero, res, niter, iters] = newton(f, df, -pi/2, 1e-10, 100)
[zero, res, niter, iters] = newton(f, df, pi, 1e-10, 100)

%%
% (c)
%
% Use fixpoint method for the negative root.

phi = @(x) x-2*f(x)/df(x);
[p, res, niter, iters] = fixpoint(phi, -pi/2, 1e-10, 100)