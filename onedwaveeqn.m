%Clear all variables
clear 
clc

%Initialise variables
L = 1;         % x in (0,L)
T = 6;       % t in (0,T)
c = 1;         % wave speed

%boolean variables that determine the nature of the simulation
triangularIC=true;
Neumann=true;
Robin=false;
damping = true;  % Set to true for damped wave equation, false for undamped

% Damping parameters (only used if damped)
alpha = 1; % spring constant
beta = 1.7;  % damping coefficient


%make the mesh
N = 200;       % cut space into N sections
M = 12000;      % cut time into M sections
dx = L/N;      % spatial grid spacing
dt = T/M;      % time step

s = c*dt/dx;   % Courant number

% Position of nodes
x = linspace(0, L, N+1);

% Initial Conditions
if triangularIC
    u0=0*x;
    a = 1/5;  % base width of the triangular displacement
    idx = find(x >= (L-a)/2 & x <= (L+a)/2); %triangular displacement
    u0(idx) = 1 - abs((x(idx) - L/2)/(a/2));  % triangular displacement

else
    u0 = sin(2*pi*x);   % initial displacement (Sin wave)

end


v0 = zeros(1,N+1);  % initial velocity

% Initialize solution arrays
u = zeros(M+1, N+1);
u(1,:) = u0;


if damping
    % Damped wave equation
    u(2,:) = u0 + dt*v0 + 0.5*s^2*(circshift(u0,-1) - 2*u0 + circshift(u0,1)) - 0.5*beta*dt*v0;
    
    % Explicit Finite Difference Scheme
    for n = 2:M
        for i = 2:N
            u(n+1,i) = (2*u(n,i) - u(n-1,i) + s^2*(u(n,i+1) - 2*u(n,i) + u(n,i-1)) + 0.5*beta*dt*u(n-1,i)) / (1 + 0.5*beta*dt);
        end
        
        % Boundary Conditions
        u(n+1,1) = 0;    % left boundary
        if Neumann
            u(n+1,N+1) = u(n+1,N);         % right boundary (Neumann)
        elseif Robin
            % Right boundary (damped)
            return
            u(n+1,N+1) = (2*u(n,N+1) - u(n-1,N+1) + s^2*(2*u(n,N) - 2*u(n,N+1)) - 2*beta*dx/dt*u(n,N+1)) / (1 + 2*beta*dx/dt + 2*alpha*dx^2/dt^2);
        else    
            u(n+1,N+1) = 0;  % right boundary (Dirichlet)
            
        end
    end
else
    u(2,:) = u0 + dt*v0 + 0.5*s^2*(circshift(u0,-1) - 2*u0 + circshift(u0,1));
    % Explicit Finite Difference Scheme
    for n = 2:M
        for i = 2:N
            u(n+1,i) = 2*u(n,i) - u(n-1,i) + s^2*(u(n,i+1) - 2*u(n,i) + u(n,i-1));
        end
        
        % Boundary Conditions
        u(n+1,1) = 0;    % left boundary
        if Neumann
            u(n+1,N+1) = u(n+1,N);         % right boundary (Neumann)
        elseif Robin
            % Right boundary (damped)
            u(n+1,N+1) = (2*u(n,N+1) - u(n-1,N+1) + s^2*(2*u(n,N) - 2*u(n,N+1)) - 2*beta*dx/dt*u(n,N+1)) / (1 + 2*beta*dx/dt + 2*alpha*dx^2/dt^2);
        else    
            u(n+1,N+1) = 0;  % right boundary (Dirichlet)
            
        end
    
    end
end    

figure(1)
betterplots
plot(x,u(1,:))
xlabel('x');
ylabel('u(x,t)');
title('Initial Displacement')
grid on
%%

if triangularIC==false
    
    % Animation
    figure(2);
    for n = 1:10:M+1
        plot(x, u(n,:), 'LineWidth', 2);
        xlabel('x');
        ylabel('u(x,t)');
        title(sprintf('1D Wave Equation (t = %.3f)', (n-1)*dt));
        ylim([-1.2 1.2]);
        grid on;
        drawnow;
        
        if n < M+1
            clf;
        end
    end
end
  
%% Animation Showing for triangle peak that also shows freq domain

if triangularIC==false
    return % do not run if not using the triangle initial condition
end 




% Animation and measurements
figure(4);
set(gcf,'units','points','position',[50,50,2000,1000])
for n = 1:10:floor(M/4)
    subplot(1,2,1);
    plot(x, u(n,:), 'LineWidth', 2);
    xlabel('x');
    ylabel('u(x,t)');
    title(sprintf('1D Wave Equation (t = %.3f)', (n-1)*dt));
    ylim([-1.2 1.2]);
    grid on;
    
    % Compute Fourier series coefficients
    U = fft(u(n,:)) / (N+1);
    U = U(1:floor((N+1)/2)+1);
    U(2:end-1) = 2*U(2:end-1);
    freq = (0:floor((N+1)/2)) / L;
    
    subplot(1,2,2);
    bar(freq, abs(U), 'LineWidth', 1.5);
    xlabel('Frequency');
    ylabel('Amplitude');
    title('Fourier Series Coefficients');
    grid on
    ylim([0 0.3])
    
 
    
    drawnow;
    
    if n < M+1
        clf;
    end
end
