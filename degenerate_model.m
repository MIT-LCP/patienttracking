function [tm,x] = degenerate_model(varargin)
%   [tm,x] = LinearModel(N,A,x0)
%
% dx=Ax
%
%
% %Degenerate mode
% A=[lamba b;0 lambda];
inputs={'N','A','x0'};
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

DET=det(A);
[v,sig]=eig(A);
tau=trace(sig);

for n=2:N
    %Apply Runge-Kutta integration method
    k1=( A*x(n-1,:)' )*dt;
    k2=( A*(x(n-1,:)'+0.5*k1) )*dt;
    k3=( A*(x(n-1,:)'+0.5*k2) )*dt;
    k4=( A*(x(n-1,:)'+k3) )*dt;
    x(n,:)=x(n-1,:)+(k1+2*k2+2*k3+k4)'./6 +[1 0];
end



