function f = bezierDCurve(n)
% Returns the derivative of a Bezier curve with degree n [(n+1) ctrl points]
f = {};
for i = 1:n+1
    f{i} = bezierD(i-1,n); 
end
end
