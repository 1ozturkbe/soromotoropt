% Adding multi_indices to path
addpath('../multi_indices');

% Problem definition in cylindrical coordinates (order P periodic bezier f'ns) 
P = 5; nKnots = 12; ptsPerKnot = 12;
[A,Ad,thetas] = influenceMatrixPeriodic(P, nKnots, ptsPerKnot);

% Creating funky function
rTest = 0.5 - sin(thetas).*(thetas-pi);
rDTest = cos(thetas).*(thetas-pi) + sin(thetas);
polarplot(thetas,rTest);

% Calculating least squares fit
sol = lsqlin(A,rTest);

% Plotting results
figure(6)
polarplot(thetas,rTest)
hold on 
polarplot(thetas,A*sol);
title('Actual function and b-spline approximation')

figure(7)
polarplot(thetas, rDTest);
hold on
polarplot(thetas, Ad*sol);
title('Actual function derivatives and b-spline approximated derivatives')

        

