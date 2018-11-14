% Adding multi_indices to path
addpath('../multi_indices')

% Example of 1-D estimation of a curve using Bezier curves
nx = 100;
x = linspace(-pi,pi,nx);
xBez = linspace(0,1,nx);
y = x.^2.*sin(x);

% Using pth order Bezier curve
p = 4;
f = bezierCurve(p);
beziervals = zeros(nx,p+1);
for i = 1:nx
    for j = 1:p+1
        beziervals(i,j) = f{j}(xBez(i));
    end
end

% Solving least squares problem
P = lsqlin(beziervals,y);
yest = beziervals*P;

% Plotting
figure(1);
plot(x,y);
hold on
plot(x,yest);