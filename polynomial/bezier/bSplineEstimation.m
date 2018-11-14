% Adding multi_indices to path
addpath('../multi_indices');

% Problem definition in cylindrical coordinates (order P bezier f'ns) 
P = 3;
nctrl = 10; nknots = nctrl; evalDeg = 1;
thetaTest = linspace(0,2*pi,100);
knots = linspace(0,2*pi,nknots);
% Creating funky function
rTest = -0.5 + sin(thetaTest).*(thetaTest-pi);
polarplot(thetaTest,rTest);

% nknots bezier f'ns, with plotting
thetaTest = linspace(0,2*pi,evalDeg*nknots);
rTest = -0.5 + sin(thetaTest).*(thetaTest-pi);
bFuncs = bezierCurve(nknots-P);
bDFuncs = bezierDCurve(nknots-P);
thetaSpacing = 2*pi/nknots;

for i = 1:size(bFuncs,2)
    thetaRange = linspace(thetaSpacing*(i-1), thetaSpacing*(i+2),20);
    valRange = linspace(0,1,20);
    fvals = bFuncs{i}(valRange);
    fDvals = bDFuncs{i}(valRange);
    figure(2);
    plot(thetaRange,fvals);
    hold on
    figure(3);
    plot(thetaRange,fDvals);
    hold on
end
figure(2);
title('Bernstein polynomials');
figure(3);
title('1st derivative of Bernstein polynomials');

% Pre-evaluating Bezier functions. 
beziervals = zeros(nknots*evalDeg, size(bFuncs,2));
bezierDs = zeros(nknots*evalDeg, size(bFuncs,2));
for i = 1:size(bFuncs,2)
    thetaRange = linspace(thetaSpacing*(i-1), thetaSpacing*(i+2),evalDeg*3);
    valRange = linspace(0,1,evalDeg*3);
    beziervals((i-1)*evalDeg+1:(i+2)*evalDeg, i) = bFuncs{i}(valRange);
    bezierDs((i-1)*evalDeg+1:(i+2)*evalDeg, i) = bDFuncs{i}(valRange);
end

% Calculating least squares fit
sol = lsqlin(beziervals,rTest);

% Recalculating interpolations and plotting
% ptsPerFn = 10;
% theta = linspace(0,2*pi,ptsPerFn*size(bFuncs,2));
% r = -0.5 + sin(theta).*(theta-pi);
% beziervals = zeros(ptsPerFn*size(bFuncs,2), size(bFuncs,2));
% bezierDs = zeros(ptsPerFn*size(bFuncs,2), size(bFuncs,2));
% for i = 1:size(bFuncs,2)
%     thetaRange = linspace(thetaSpacing*(i-1), thetaSpacing*(i+2),ptsPerFn);
%     valRange = linspace(0,1,evalDeg*3);
%     beziervals((i-1)*evalDeg+1:(i+2)*evalDeg, i) = bFuncs{i}(valRange);end
% figure(5);
% polarplot(theta,r);

% Plotting results
figure(6)
polarplot(thetaTest,rTest)
hold on 
polarplot(thetaTest,beziervals*sol);
title('Actual function and b-spline approximation')
%polarplot(theta,1./bezierDs*sol);

        

