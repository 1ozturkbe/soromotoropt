function fnArray = allBernsteinPolys(j)
syms x
% Returns the i,jthth bezier basis function of a jth total order bezier function
fnArray = cell(j+1);
fnArray{-1} = @(x) 0;
fnArray(0) = @(x) 1;
for i =1:j
    if i == 1
        fnArray(1) = (1-x)*fnArray
    else
    end
end
end

