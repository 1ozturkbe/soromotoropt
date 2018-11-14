% Adding multi_indices to path
addpath('../multi_indices')

% 2D problem definition (estimated by product of Bezier curves on x,y);
nx = 50; ny = 50; r = 1;
x = linspace(-pi,pi,nx);
xBez = linspace(0,1,nx);
y = linspace(-r,r,ny);
yBez = linspace(0,1,ny);
[X,Y] = meshgrid(x,y);
C = @(x,y) cos(x).*y.^2;
Ceval = C(X,Y);
surf(X,Y,Ceval);

% Order for bezier approx
Px = 10; Py = 10;

% Getting the proper multi-indices (args are degree, and total order)
tomil = TotalOrderMultiIndexLattice(2, max(Px,Py));
tomil.init();
Mtotalorder = tomil.get_midx_matrix();
invalidMIs = [];

% Truncating elements
for i = 1:size(Mtotalorder,1)
    if Mtotalorder(i,1) > Px
        invalidMIs = [invalidMIs, i];
    elseif Mtotalorder(i,2) > Py
        invalidMIs = [invalidMIs, i];
    end
end
Mtotalorder(invalidMIs,:) = [];

% Setting up 2D Bezier curves
Cx = {};
Cy = {};

fxint = bezierCurve(Px);
fyint = bezierCurve(Py);
beziervalsx = zeros(nx,Px+1);
beziervalsy = zeros(ny,Py+1);
for ix = 1:nx
   for px = 1:size(fxint,2)
      beziervalsx(ix,px) = fxint{px}(xBez(ix));
   end
end
for iy = 1:ny
    for py = 1:size(fyint,2)
        beziervalsy(iy,py) = fyint{py}(yBez(iy));
    end
end   

% Regularizing and flattening problem
Cevalflat = reshape(Ceval, [size(Ceval,1)*size(Ceval,2),1]);
Amat = zeros(size(Cevalflat,1),size(Mtotalorder,1));
for i = 1:size(Amat,1)
    for j = 1:size(Amat,2)
        ix = ceil(i/nx); iy = mod(i,nx)+1;
        ox = Mtotalorder(j,1); oy = Mtotalorder(j,2);
        Amat(i,j) = sum(beziervalsx(ix,ox+1).*beziervalsy(iy,oy+1));
    end
end

% Solving the optimization problem
sol = lsqlin(Amat,Cevalflat);

% Plotting results
figure(3);
subplot(1,3,1);
surf(x,y,Ceval)
subplot(1,3,2);
surf(x,y,reshape(Amat*sol,[nx,ny]))
subplot(1,3,3)
surf(x,y,abs(Ceval-reshape(Amat*sol,[nx,ny])))


