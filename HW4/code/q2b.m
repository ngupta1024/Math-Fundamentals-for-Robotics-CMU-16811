clear;

x0 = [1, -1];
step = main(x0);

function step = main(x)
    [dx, dy] = gradients(x);
    step = 0;
    path = [x];
    while abs(dx) > 10^-5 || abs(dy) > 10^-5
        t = cal_t(x, dx, dy);
        x = x - t*[dx, dy];
        path = [path; x];
        [dx, dy] = gradients(x);
        step = step + 1;
    end
    plot_path(path);
end

function [dx, dy] = gradients(x)
    dx = 3*x(1)^2-4*x(1);
    dy = 3*x(2)^2+6*x(2);
end

function t = cal_t(x, dx, dy)
    syms tt;
    dg = -3*dx*(x(1)-tt*dx)^2-3*dy*(x(2)-tt*dy)^2+4*dx*(x(1)-tt*dx)-6*dy*(x(2)-tt*dy);
    result = solve(dg, tt);
    if result(1) < 0 && result(2) > 0
        t = double(result(2));
    elseif result(1) > 0 && result(2) < 0
        t = double(result(1));
    elseif -3*dx^3-3*dy^3 < 0
        t = double(result(1));
    else
        t = double(result(2));
    end
end

function plot_path(path)
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
    plot(dX, dY, 'b*');
    hold on;
    contour(X, Y, Z, c, 'ShowText','on');
    hold on;
    [xt,yt] = meshgrid(-3:0.5:3,-4:0.5:3);
    zt = -(xt.^3 + yt.^3 - 2*xt.^2 + 3*yt.^2 - 8);
    [gx, gy] = gradient(zt);
    quiver(xt,yt,gx,gy);
    
    drawArrow = @(p1,p2,varargin) quiver( p1(1),p1(2),p2(1)-p1(1),p2(2)-p1(2),0, varargin{:} );
    for p = 2:length(path)
        p1 = path(p-1, :);
        p2 = path(p, :);
        drawArrow(p1,p2,'linewidth',2,'color','r'); 
        hold on;
    end
    hold off;
end