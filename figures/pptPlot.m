nx = 12;
nt = 4;

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

% Mapping Chebyshev polynomials
close all
d = 0:3;
n = 100;
theta = -1:0.01:1;
for i = d
    plot(pi*theta, 1*(chebyshevT(2*i, theta)+0),'LineWidth', 2);
    hold on
end

% Example result with Chebyshev polys
close all
polarplot(pi*theta, -0.075*(chebyshevT(6, theta)+0.5) + ...
    0.1*(chebyshevT(8, theta)+1.25) + 0.25*(chebyshevT(2, theta)+1.5),'LineWidth', 2);

% Mapping periodic sinusoids
close all
d = 2:5;
n = 100;
theta = -1:0.01:1;
for i = d
    polarplot(pi*theta, sin(pi*theta*i/2));
    hold on
end
% Example result with sinusoids
close all
polarplot(pi*theta, (2+ sin(pi*theta*3/2) + sin(pi*theta*5/2))/4);
