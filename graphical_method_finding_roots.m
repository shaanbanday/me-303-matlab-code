% set x as a discretized data set with intervals of 0.01 from -10 to 10
x = [-10:0.01:10];

f = sin(x); % f(x) = sin (x), represented as an array of discretized values
g = cos(pi*x); % g(x) = cos(π*x), represented same way as f(x)

y = f - g; % y = f(x) - g(x)

%plot figures
subplot (2, 1, 1)
plot(x, y, LineWidth= 2, LineStyle="-");
hold on
plot (x, zeros(size(x)), '--k')
legend("sin(x) - cos(\pix)")
ylabel("y")
xlabel("x")

subplot (2,1,2)
plot(x,f, LineWidth= 2)
hold on
plot (x,g,LineWidth=2, LineStyle= "-")
legend("f(x) = sin(x)", "g(x) = cos(\pix)")
ylabel("y")
xlabel("x")