function [ phase ] = estimate_phase_portrait(x,show)
%[ phase ] = estimate_phase_portrait(x,show)
%
% Estimates phase portrait of two variables (x is Nx2).
% Assumes time sampling of x is constant.
%


%Estimate phase space -maybe estimate using Runge-Kutta?
phase(:,1)=[diff(x(:,1)) ;NaN];
phase(:,2)=[diff(x(:,2));NaN];

if(show)
    figure
    quiver(x(:,1),x(:,2),phase(:,1),phase(:,2))
    xlabel('x1');ylabel('x2')
    grid on
end
end

