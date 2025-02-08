clear; clc;

% inputs
A = [4, -1, 0;
    -1, 4, -1;
    0, -1, 4];
b = [3;5;6];

x0 = zeros(size(b));        % initial guess (prev guess)
x = zeros(size(b));         % store updated guess
TOL = 1e-5;
N = 100;                    % max iterations
n = length(b);
alpha = 0.75;                  % relaxation method coeff. (<1 = under-relaxation)

% G-S
%S1
k = 1;

%S2
while k<=N

    %S3
    for i=1:n
        sum1 = 0;
        sum2 = 0;

        for j=1:n
            if j < i
                sum1 = sum1 + A(i,j)*x(j);
            elseif j > i
                sum2 = sum2 + A(i,j)*x0(j);
            end
        end

        % update x with relaxation method
        x(i) = alpha * ((-sum1 - sum2 + b(i)) / A(i,i)) + (1-alpha)*x0(i);     
    end

    %S4
    if norm(x-x0, 2) < TOL
        disp(x);
        fprintf('Converged after %d iterations', k);
        break;
    end

    %S5
    k = k+1;                        % update k

    %S6
    for i=1:n
        x0(i) = x(i);               % set prev guess = updated guess
    end
end

if k>N
    fprintf('Did not converge');
end




