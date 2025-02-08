clear; clc;

a = 0;           % starting point
b = 5;           % endpoint
dx = 0.01;
x = a:dx:b;    % discretize domain

% define function
f = @(g) (2*g.^3 - 12);      % function f, given input g
df = @(g) (6*g.^2);  % derivative of f

y = f(x);                       % compute y = f(x)

y_zeros = zeros(size(x));

% could find the root graphically
plot(x,y, LineWidth=3);
hold on;
plot(x,y_zeros,'k', LineWidth=3);        % plot a line at y=0
grid on;
ylabel('y');
xlabel('x');

% NEWTON RAPHSON METHOD
p0 = 5;   % initial (prev) guess
p = 0;    % next guess
TOL = 1e-6;
N = 100;  % max iterations

% step1
i =1;

% step 2
while i<=N
    %step 3
    p = p0 - f(p0)/df(p0);

    % step 4
    if abs(p-p0) < TOL
        disp(p);
        fprintf('Converged in %d iterations', i);
        break;
    end
    
    % step 5
    i=i+1;

    % step 6
    p0 = p;

end


if i>N
    fprintf('Did not converge');
end
