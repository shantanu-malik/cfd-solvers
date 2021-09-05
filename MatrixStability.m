% constructing matrix A
c = zeros(100);
c(1,[1,2]) = [1, -1];
c(end,[end-1,end]) = [1, -1];
for j = 2:99
    c(j,j-1) = 1/2;
    c(j,j+1) = -1/2;
end

% contructing Euler stability region (circle)
x = zeros(1,361); y=x;
for t = 0:360
    x(t+1) = cos(deg2rad(t)) - 1;
    y(t+1) = sin(deg2rad(t));
end

plot(x,y)
hold on
plot(eig(c))
xlim([-3 1])
ylim([-2 2])
daspect([1 1 1])