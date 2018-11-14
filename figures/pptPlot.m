nx = 24;
nt = 6;

r = zeros(nx+1,nt);
theta = zeros(nx+1,nt);

for i = 1:nx+1
    for j = 1:nt
        theta(i,j) = 2*pi/nx*(i+0.2*j);
        if mod(i,3) > 0
            r(i,j) = 1/(nt/j)^1.2;
        else
            r(i,j) = 1/(0.5 + 0.5*nt/j)^1.2;
        end
    end
end

polarplot(theta',r', '-')
hold on
polarplot(theta, r, '.-')