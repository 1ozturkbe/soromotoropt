function f = bezierCurve(n)
% Returns a Bezier curve with degree n [(n+1) ctrl points]
f = {};
for i = 1:n+1
    f{i} = bezier(i-1,n); 
end
end
