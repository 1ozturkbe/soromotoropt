function d = bezierD(i,n)
% Returns the derivative of the ith function of an n-term Bezier curve
syms t
if i == 0
    a = @(t) 0;
    b = bezier(i,n-1);
elseif i == n
    a = bezier(i-1,n-1);
    b = @(t) 0;
else
    a = bezier(i-1,n-1);
    b = bezier(i,n-1);
end
    d = @(t) n*(a(t)-b(t));
end

