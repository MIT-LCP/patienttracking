function [tm,x] = linearmodel(varargin)
%   [tm,x] = LinearModel(N,A,x0,show)
%
% dx=Ax
%
%
% %Degenerate mode
% A=[
inputs={'N','A','x0','show'};
N=[];
x0=[7 20]; %Default initial conditions
A=[1 0; 0 1]; %Default State space matrix

for n=1:nargin
    if(~isempty(varargin{n}))
        eval([inputs{n} '=varargin{n};'])
    end
end

dt=1/10; %Sampling interval in hours
tm=linspace(0,dt*(N-1),N);
x=zeros(N,2);
x(1,:)=x0(:)';


for n=2:N
    %Apply Runge-Kutta integration method
    k1=( A*x(n-1,:)' )*dt;
    k2=( A*(x(n-1,:)'+0.5*k1) )*dt;
    k3=( A*(x(n-1,:)'+0.5*k2) )*dt;
    k4=( A*(x(n-1,:)'+k3) )*dt;
    x(n,:)=x(n-1,:)+(k1+2*k2+2*k3+k4)'./6 +[1 0];
end

if(show)
    rge=min(x(:,1)):max(x(:,1));
    %Get eigenvectors and eigenvalues and plot the in phase space
    [v,s]=eig(A);
    s1=v(2,1)/v(1,1);
    s2=v(2,2)/v(1,2);
    plot(rge,rge*s1,'k-');hold on;
    plot(rge,rge*s2,'k-')
    
end

