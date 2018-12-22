clear;

x = -3:0.05:3;
y = -4:0.05:3;
[X,Y] = meshgrid(x,y);
Z = X.^3 + Y.^3 - 2*X.^2 + 3*Y.^2 - 8;

dx = [0, 4/3];
dy = [0, -2];
[dX, dY] = meshgrid(dx,dy);
c = dX.^3 + dY.^3 - 2*dX.^2 + 3*dY.^2 - 8;
c = sort(c(:)');

figure();
contour(X, Y, Z, 'ShowText','on');
hold on;
contour(X, Y, Z, c, 'ShowText','on');
hold on;
plot(dX, dY, 'b*');
hold on;
[xt,yt] = meshgrid(-3:0.5:3,-4:0.5:3);
zt = -(xt.^3 + yt.^3 - 2*xt.^2 + 3*yt.^2 - 8);
[gx, gy] = gradient(zt);
quiver(xt,yt,gx,gy);
hold off;
