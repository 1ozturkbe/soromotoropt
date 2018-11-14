function [A,Ad,thetas] = influenceMatrix(P,nKnots, ptsPerKnot)
% Return the influence matrix of an order P b-spline function, 
% with nKnots and ptsPerKnot
bFuncs = bezierCurve(nKnots-P);
bDFuncs = bezierDCurve(nKnots-P);
thetaSpacing = 2*pi/nKnots;
A = zeros((nKnots-P+1)*ptsPerKnot, size(bFuncs,2));
Ad = zeros((nKnots-P+1)*ptsPerKnot, size(bFuncs,2));
thetas = linspace(0,2*pi,nKnots*ptsPerKnot);
for i = 1:size(bFuncs,2)
    thetaRange = linspace(thetaSpacing*(i-1), thetaSpacing*(i+P-1),P*ptsPerKnot);
    valRange = linspace(0,1,P*ptsPerKnot);
    fvals = bFuncs{i}(valRange);
    fDvals = bDFuncs{i}(valRange);
    A((i-1)*ptsPerKnot+1:(i+P-1)*ptsPerKnot,i) = fvals;
    Ad((i-1)*ptsPerKnot+1:(i+P-1)*ptsPerKnot,i) = fDvals;
end