clear all
clc;

x = linspace(0, 2*pi, 100); %100 equally spaced points from 0 to 2 pi

y = x; %100 equally spaced points from 0 to 2 pi

[X, Y] = meshgrid(x, y);

Z = sin(X) .* cos(Y);

s = surf(X, Y, Z);

s.EdgeColor = "none";

colorbar;

colormap = "jet"

title("Tutorial 1 3D Plot")