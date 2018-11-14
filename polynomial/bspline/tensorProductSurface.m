% Eg. Using Berstein poly approx. in 2D to create a random surface.
syms x y
b3 = bernstein(sin(2*pi*x), 3, x);
e3 = bernstein(exp(y), 3, y);
% Plotting approximations
figure(1);
fplot(sin(2*pi*x),[0,1])
hold on
fplot(b3,[0,1])
fplot(e3,[0,1])

% Plotting tensor product of the two functions
figure(2);
xgrid = linspace(0,1,100);
ygrid = xgrid;
[X,Y] = meshgrid(x,y);
fsurf(e3*b3,[0,1,0,1])