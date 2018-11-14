function f = bezier(i,n)
% Returns the ith function of an n-term Bezier curve
f = @(t) nchoosek(n,i) * t.^i .* (1-t).^(n-i);
end

