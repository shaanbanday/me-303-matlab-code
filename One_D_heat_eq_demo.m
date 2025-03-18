%Clear all variables
clear 
clc
%Initialise variables

L= 1;         % x in (0,L)
T= 0.5;       % t in (0,T)

k=1;    % conductivity

%make the mesh
N=30;   % cut space into N sections
M=5000; % cut time  into M sections
dx=L/N; dt=T/M; % grid spacing

F=k*dt/dx^2;
 
% Position of nodes
for i=1:N+1
x(i)=(i-1)*dx; %find positions of x nodes
end
 
% Initial Condition
for i=1:N+1
T0(i)=sin(2*pi*x(i)); %Set initial temp dist
end
 
% Explict Scheme for Partial Difference Equation
for j=1:M % time corrdinate = j/M
    
    for i=2:N % space corrdinate = i/N
        T1(i)=T0(i)+F*(T0(i+1)-2*T0(i)+T0(i-1)); %internal nodes
    end

    %Boundary conditions
    T1(1)=1; % DBC left
    
%   T1(N+1)=0; % DBC right: a constant BC
    T1(N+1)= 5; % DBC right: a time-varying one

    T0=T1;
    Temp(j,:)=T1;
end
 
%% plot
figure(1)
[X,Y] = meshgrid(0:dx:L,dt:dt:T); 
mesh(X,Y,Temp); 
  shading interp
colormap('jet')
xlabel('x'); ylabel('t'); zlabel('T(x,t)'); colorbar
title("Temperature distribution over time")

%% plot
figure(2)
plot(Temp(1:1000:end, :)')
xlabel('x'); 
ylabel('T(x,t)');
grid on
title("Temperature Distribution every 1000 timesteps")
