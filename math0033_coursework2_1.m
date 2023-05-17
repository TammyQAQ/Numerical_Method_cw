%% Yitian (Tammy) Wang
% MATH0033 Numerical Methods Computational homework 1
clear all, close all, clc, format long, format compact
%%
% *Exercise 1*
%
% (a)
%
% Strictly diagonal dominant when sum of a_ij where i!=j is smaller than
% a_ii, for this to hold for A_epsi, 1>2*epsilon+2*epsilon^2 as espilon is between 0
% and 1. Solve the inequality to get the bound for epsilon,
% (-sqrt(3)-1)/2<epsilon<(sqrt(3)-1)/2. As epsilon is between 0 and 1, we have
% 0<epsilon<(sqrt(3)-1)/2
% 
%%
% (b)

n=5;
epsi=0.3;
[A,b]=matrix(n,epsi);
x0=zeros(n,1);
tol=1e-10;
nmax=1e3;
[Jx,Jiter]=itermeth(A,b,x0,nmax,tol,'J')
[Gx,Giter]=itermeth(A,b,x0,nmax,tol,'G')

% Jacobi converges in 50 iterations and Gauss_Seidel converges in 14
% iterations.

%%
% (c)

% set epsi to 0, 0.01, 0.02,...,1.
n=5;
N=101;
epsi=linspace(0,1,N);
rho_J=zeros(N,1);
rho_GS=zeros(N,1);
for i=1:N
    A=matrix(5,epsi(i));
    D=diag(diag(A));
    L=tril(A)-D;
    U=triu(A)-D;
    B_J=-D^(-1)*(L+U);
    B_GS=-(L+D)^(-1)*(U);
    rho_J(i)=max(abs(eig(B_J)));
    rho_GS(i)=max(abs(eig(B_GS)));
    disp(sprintf('for epsilon=%g',i/100));
    disp(sprintf('Spectral radius for Jacobi is %g',rho_J(i)));
    disp(sprintf('Spectral radius for Gauss-Siedel is %g',rho_GS(i)));
end
plot(epsi,rho_J,epsi,rho_GS);
hold on;
xlabel('epsilon');
ylabel('Spectrual radius');
legend('Jacobi','Gauss-Seidel');

%%
% Spectral radius of Jacobi converges for epsilon<0.45 and of Gauss
% converges for epsilon<0.71.
% Comparing the result to (a), Jacobi agreed with the bound
% 0<epsilon<(sqrt(3)-1)/2, where (sqrt(3)-1)/2 approximately equals to
% 0.37. However, Gauss-Siedel's spectral radius is much larger, thus it
% converges quicker. This agrees with (b) as well.

%%
% For n=5, epsilon=0.5, the iteration matrix will be greater than 1 for
% Jacobi but less than 1 for Gauss, thus Gauss will converge where as
% Jacobi will not. 

n=5;
[A,b]=matrix(5,0.5);
x0=zeros(5,1);
[x,niter_J]=itermeth(A,b,x0,1e3,1e-10,'J')
[x,niter_G]=itermeth(A,b,x0,1e3,1e-10,'G')

%%
% From the above result, we see that Jacobi reached the maximum iteration
% number without converging. Gauss converged in 27 iterations. The outcome
% consist with my recommendation.









