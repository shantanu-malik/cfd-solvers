function solver2()
%% 2nd order central scheme (interior) with 1st order boundary schemes
% This function solves the convection equation with c = 1 for a sinusoidal
% function using finite difference schemes.
% Spacial scheme: 2nd order central (interior) with 1st order boundary schemes
% Temporal scheme: Forward Euler scheme
% Stability: Unstable for all dt

%% Initialization
T = 20*pi;   % simulation time
dx = 2*pi/100;  % spacial grid size
dt = 0.01*dx; % temporal grid size
x = 0:dx:2*pi;  % spacial domain
t = 0:dt:T; % temporal domain

u = zeros(length(t),length(x));
error = zeros(1,length(t));
u_max = error;
u(1,:) = sin(x);    % initial condition
u_max(1) = max(u(1,:));

%% Solver
for n = 2:size(u,1)
    for j = 1:size(u,2)
        if j == 1
            u_x = (u(n-1,j+1) - u(n-1,j))/(dx);
        elseif j == size(u,2)
            u_x = (u(n-1,j) - u(n-1,j-1))/(dx);
        else
            u_x = (u(n-1,j+1) - u(n-1,j-1))/(2*dx);
        end
        u_t = - u_x;
        u(n,j) = u(n-1,j) + dt*u_t;
    end
    % Error and u_max calculation
    u_max(n) = max(abs(u(n,:)));
    error(n) = (rms(u(n,:) - sin(x-t(n))))*100/max(sin(x-t(n)));
end

%% Plotting the results
% Solution
figure
tic
for n = 1:max(1,floor(0.1/dt)):size(u,1)
    pause((n-1)*dt - toc)
    subplot(2,2,[1,2])
    plot(x,sin(x-t(n)),'--k',x,u(n,:),'r');grid on
    xlim([0,2*pi]);ylim([-15,15])
    title(['Time:',num2str((n-1)*dt),' | Red: Numerical | Black: Analytical']);
    xlabel('X');ylabel('u')
end

% Error
subplot(2,2,3)
plot(t,error)
grid on;xlim([0,T])
title('Error_{RMS}');xlabel('t');ylabel('Error %')
subplot(2,2,4)
plot(t,u_max)
grid on;xlim([0,T])
title('u_{max}');xlabel('t');ylabel('u_{max}')
end