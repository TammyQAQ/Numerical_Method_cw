%% Yitian (Tammy) Wang
% MATH0033 Numerical Methods Computational homework 1
clear all, close all, clc, format long, format compact
%%
% *Exercise 1*
%
% (a)

% Define the function and exact solution
f = @(t,y)1-y^2;
y_ex = @(t)tanh(t);
dfdy = @(t,y)-2*y;
% Set the initial condition
y0 = 0;
% Set the time interval and the number of time steps
tmax = 20;
% Iterate through N and plot graph
[t_f19, u_f19] = feuler(f, [0,tmax], y0, 19)
[t_f21, u_f21] = feuler(f, [0,tmax], y0, 21)
[t_f40, u_f40] = feuler(f, [0,tmax], y0, 40)
[t_f80, u_f80] = feuler(f, [0,tmax], y0,80)
[t_f160, u_f160] = feuler(f, [0,tmax], y0, 160)
t = tmax*(0:0.01:1);      % t values for plotting exact solution
y = y_ex(t);
plot(t, y, 'k', t_f19, u_f19, t_f21, u_f21', t_f40, u_f40, t_f80, u_f80, t_f160, u_f160)
legend('exact','N=19','N=21','N=40','N=80','N=160')

%%
% Compute the errors e_fe
e_fe = zeros(5,1);
Nvec = [19, 21, 40, 80, 160];
for i = 1:5
    [t_fi, u_fi] = feuler(f, [0,tmax], y0, Nvec(i));
    y = y_ex((0:length(u_fi)-1)*20/Nvec(i));
    e_fe(i) = max(abs(y-u_fi))
end

%%
% Compute the errors e_feend
e_feend = zeros(5,1);
Nvec = [19, 21, 40, 80, 160];
ytmax = y_ex(tmax);
for i = 1:5
    N = Nvec(i);
    [t_fi, u_fi] = feuler(f, [0,tmax], y0, Nvec(i));
    e_feend(i) = abs(ytmax - u_fi(end))
end

%%
% (b)

% Define the function and exact solution
f = @(t,y)1-y^2;
y_ex = @(t)tanh(t);
dfdy = @(t,y)-2*y;
% Set the initial condition
y0 = 0;
% Set the time interval and the number of time steps
tmax = 20;
% Iterate through N and plot graph
[t_h19, u_h19] = heun(f, [0,tmax], y0, 19)
[t_h21, u_h21] = heun(f, [0,tmax], y0, 21)
[t_h40, u_h40] = heun(f, [0,tmax], y0, 40)
[t_h80, u_h80] = heun(f, [0,tmax], y0, 80)
[t_h160, u_h160] = heun(f, [0,tmax], y0, 160)
t = tmax*(0:0.01:1);     % t values for plotting exact solution
y = y_ex(t);
plot(t, y, 'k', t_h19, u_h19, t_h21, u_h21', t_h40, u_h40, t_h80, u_h80, t_h160, u_h160)
legend('exact','N=19','N=21','N=40','N=80','N=160')

%%
% Compute the errors e_fe
e_hune = zeros(5,1);
Nvec = [19, 21, 40, 80, 160];
for i = 1:5
    [t_hi, u_hi] = heun(f, [0,tmax], y0, Nvec(i));
    y = y_ex((0:length(u_hi)-1)*20/Nvec(i));
    e_heun(i) = max(abs(y-u_hi))
end

%%
% Compute the errors e_feend
e_heunend = zeros(5,1);
Nvec = [19, 21, 40, 80, 160];
ytmax = y_ex(tmax);
for i = 1:5
    N = Nvec(i);
    [t_hi, u_hi] = heun(f, [0,tmax], y0, Nvec(i));
    e_heunend(i) = abs(ytmax - u_hi(end))
end

%%
% (c)
% Compare e_fe and e_hune for N
Nvec = [19, 21, 40, 80, 160];
loglog(Nvec, e_fe, Nvec, e_heun)
legend('error fe', 'error heun')
xlabel('N')
ylabel('max error')

% Theory suggest that Hune's method is second order convergent where as
% Foward Euler is first order. These agree with the graph that is
% obtained.

%%
% (d)

% No, not every approximation obtained by Foward Euler and Hune's reproduce
% the same asymptotic behaviour. We can see this from garph produced in
% part a and b that when N = 19, both methods do not present asymptotic
% behaviour.
%
% Both methods approximates well for y(20), for large N, error could be
% neglected.
%
% Using h = N/20, h > 1 lack of stability and h < 1 preforms well, thus both
% methods are stable to be h = 1.
%
% The lack of stability at N = 19 for Forward Euler manifest into
% oscilation around the asymptote increases and divergeses away from exact
% solution. And for Hune's method, the asymptote is convengeing to around
% 0.6 for the exact solution where the asympotote is 1.


