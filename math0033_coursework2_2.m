% *Exercise 2*
%%
% (a)

rhoBJ=zeros(5,1);
rhoBGS=zeros(5,1);
Nvec=[5,10,20,40,80];
for i=1:5
    N=Nvec(i)
    h=1/N;
    A=(2/h^2)*diag(ones(N-1,1))-(1/h^2)*diag(ones(N-2,1),1)-(1/h^2)*diag(ones(N-2,1),-1);
    b=transpose(sin(pi*h*(1:N-1)));
    % construct D,L,U
    D=diag(diag(A));
    L=tril(A)-D;
    U=triu(A)-D;
    BJ=-D^(-1)*(L+U);
    BGS=-(D+L)^(-1)*U;
    rhoBJ(i)=max(abs(eig(BJ)));
    rhoBGS(i)=max(abs(eig(BGS)));
    fprintf('spectral radius of Jacobi is %g\n', rhoBJ(i));
    fprintf('spectral radius of Gauss-Siedel is %g\n', rhoBGS(i));
end

%%
% From the result, both method would be convergent, as A is a
% symmetric matrix and spetrial radius is smaller than 1 for all N.
% Comparing the two radius, for each N, Gauss_Siedel has a smaller radius,
% thus it will converge faster than Jacobi.
%
% As the system size N increases the spectral radius approaches to 1, the
% proformance of both method would decrease (slower convergance).

%%
Nvec=[5,10,20,40,80];
loglog(Nvec,1-rhoBJ)
hold
loglog(Nvec,1-rhoBGS)
legend('Jacobi','Gauss-Seidel')

%%
% log(1-rho) has an linear relationship with log(N), thus
% log(1-rho)=alpha*log(N)+C, where alpha and C are constant. To find C and
% alpha, it is enough to know the gradient of the above plot. C = e^b =
% e^(log(1-rho)-alpha*log(N)) = (1-rho)*N^(-alpha)

grad_J=(log(1-rhoBJ(5))-log(1-rhoBJ(4)))/(log(80)-log(40))
C=(1-rhoBJ(5))*80^(-grad_J)
grad_GS=(log(1-rhoBGS(5))-log(1-rhoBGS(4)))/(log(80)-log(40))
C=(1-rhoBGS(5))*80^(-grad_GS)

%%
% Using the assumption that the error in both methods is approximately
% proportional to rho(B)^k, we could express using kth iteration. Thus
% e^k=d(1-CN^alpha)^k, for some constant d. By taking log and using Laurent
% expansion, we deduce that tolerence cost when N large for both method is
% O(N^3). Which is less efficient than direct method where the cost is O(N). 

%%
% (b)
clear all, close all, clc, format long, format compact
Nvec=[5,10,20,40,80];
for i=1:5
    N=Nvec(i)
    h=1/N;
    A=(2/h^2)*diag(ones(N-1,1))-(1/h^2)*diag(ones(N-2,1),1)-(1/h^2)*diag(ones(N-2,1),-1);
    b=transpose(sin(pi*h*(1:N-1)));
    [u, niter]=itermeth(A, b, zeros(N-1, 1), 1e5, 1e-10, 'G')
    x=linspace(0, 1, N+1);
    plot(x, [0;u;0])
    hold on
end
plot(x,pi^(-2)*sin(pi*x))
legend('N=5','N=10','N=20','N=40','N=80','exact solution')

%%
% from the plot, we could see that as N increase the accuracy increace,
% i.e. the line is closer to the exact solution.

%%
clear all, close all, clc, format long, format compact
Nvec=[5,10,20,40,80];
error_vect=zeros(5,1);
for i=1:5
    N=Nvec(i);
    h=1/N;
    A=(2/h^2)*diag(ones(N-1,1))-(1/h^2)*diag(ones(N-2, 1), 1)-(1/h^2)*diag(ones(N-2,1),-1);
    b=transpose(sin(pi*h*(1:N-1)));
    [u,niter]=itermeth(A, b, zeros(N-1, 1), 1e5, 1e-10, 'G');
    u_n=[0; u(1:length(u)); 0];
    x=linspace(0, 1, N+1)';
    y=pi^(-2)*sin(pi*x);
    error_vect(i)=max(abs(u_n-y));
    fprintf('The error is %2.10f\n', error_vect(i))
end

%%
% Plot
Nvec=[5,10,20,40,80];
hvect=1./Nvec;
loglog(hvect, error_vect)

%%
p=(log(error_vect(5))-log(error_vect(4)))/(log(hvect(5))-log(hvect(4)))
C=error_vect(5)*hvect(5)^(-p)
loglog(hvect, error_vect)
hold
loglog(hvect, C*hvect.^p)
%The two lines are fitting very well, numerical results support the
%theoretical estimate.

%%
% (c)

clear all, close all, clc, format long, format compact
Nvec=[5,10,20,40,80];
error_vect=zeros(5,1);
for i=1:5
    N=Nvec(i);
    h=1/N;
    A=(2/h^2)*diag(ones(N-1,1))-(1/h^2)*diag(ones(N-2, 1), 1)-(1/h^2)*diag(ones(N-2,1),-1);
    b=ones(N-1,1);
    [u,niter]=itermeth(A,b,zeros(N-1,1),1e5,1e-10,'G')
    x=linspace(0,1,N+1);
    plot(x,[0;u;0])
    hold on
end
plot(x,(1/2)*x.*(1-x))
legend('N=5','N=10','N=20','N=40','N=80','exact solution')

%%
% Compute the error:
clear all, close all, clc, format long, format compact
Nvec=[5,10,20,40,80];
error_vect=zeros(5,1);
for i=1:5
    N=Nvec(i);
    h=1/N;
    A=(2/h^2)*diag(ones(N-1,1))-(1/h^2)*diag(ones(N-2, 1), 1)-(1/h^2)*diag(ones(N-2,1),-1);
    b=ones(N-1,1);
    [u,niter]=itermeth(A, b, zeros(N-1, 1), 1e5, 1e-10, 'G');
    u_n=[0; u(1:length(u)); 0];
    x=linspace(0, 1, N+1)';
    y=(1/2).*x.*(1-x);
    error_vect(i)=max(abs(u_n-y));
    fprintf('The error is %g\n', error_vect(i))
end

%%
% Plot:
hvect=1./Nvec;
loglog(hvect, error_vect)

%%
% From the graph, we can see that the error increases as N increase,
% however, error converges to a fixed value which is smaller than the
% tolerance. Thus method is always convergent. As the fourth derivative of
% y is 0, and C is propotional to this, thus C=0<CN^(-2).
