function solver2rk4()
%% 2nd order central scheme (interior) with 1st order boundary schemes
% This function solves the convection equation with c = 1 for a sinusoidal
% function using finite difference schemes.
% Spacial scheme: 2nd order central (interior) with 1st order boundary schemes
% Temporal scheme: RK-4 scheme
% Stability: Stable for dt <= 2.83*dx

%% Initialization
T = 20*pi;   % simulation time
dx = 2*pi/100;  % spacial grid size
dt = 2.8*dx; % temporal grid size
x = 0:dx:2*pi;  % spacial domain
t = 0:dt:T; % temporal domain

u = zeros(length(t),length(x));
k1 = zeros(1,length(x)); k2 = k1; k3 = k1; k4 = k1;
error = zeros(1,length(t));
u_max = error;
u(1,:) = sin(x);    % initial condition
u_max(1) = max(u(1,:));

%% Solver
for n = 2:size(u,1)
    U = u(n-1,:);
    for j = 1:size(u,2)
        if j == 1
            u_x = (U(j+1) - U(j))/(dx);
        elseif j == size(u,2)
            u_x = (U(j) - U(j-1))/(dx);
        else
            u_x = (U(j+1) - U(j-1))/(2*dx);
        end
        u_t = - u_x;
        k1(j) = dt*u_t;
    end
    
    U = u(n-1,:) + k1/2;
    for j = 1:size(u,2)
        if j == 1
            u_x = (U(j+1) - U(j))/(dx);
        elseif j == size(u,2)
            u_x = (U(j) - U(j-1))/(dx);
        else
            u_x = (U(j+1) - U(j-1))/(2*dx);
        end
        u_t = - u_x;
        k2(j) = dt*u_t;
    end
    
    U = u(n-1,:) + k2/2;
    for j = 1:size(u,2)
        if j == 1
            u_x = (U(j+1) - U(j))/(dx);
        elseif j == size(u,2)
            u_x = (U(j) - U(j-1))/(dx);
        else
            u_x = (U(j+1) - U(j-1))/(2*dx);
        end
        u_t = - u_x;
        k3(j) = dt*u_t;
    end
    
    U = u(n-1,:) + k3;
    for j = 1:size(u,2)
        if j == 1
            u_x = (U(j+1) - U(j))/(dx);
        elseif j == size(u,2)
            u_x = (U(j) - U(j-1))/(dx);
        else
            u_x = (U(j+1) - U(j-1))/(2*dx);
        end
        u_t = - u_x;
        k4(j) = dt*u_t;
    end
    
    u(n,:) = u(n-1,:) + (k1 + 2*k2 + 2*k3 + k4)/6;
    
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