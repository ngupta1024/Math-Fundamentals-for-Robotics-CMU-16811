clear;

interval = [-pi/2, pi/2];
plot_fx(interval);

approx = muller(-pi/2, -pi/4, 0, 2);

function y = fx(x)
    y = 0.5+sin(x);
end

function plot_fx(interval)
    x = interval(1):0.1:interval(2);
    y = fx(x);
    figure(1);
    plot(x, y);
    legend('y=0.5+sin(x)');
end

function [y] = tx(x)
    y = (pi/2-x)*cos(x) + sin(x) - 1;
end

function [y] = gx(x, roots)
    y = tx(x);
    [~, n] = size(roots);
    for i = 1:n
        y = y/(x-roots(i));
    end
end

function [a, b, c] = get_coefficient(x0, x1, x2, y0, y1, y2)
    c = y2;
    
    div_diff1 = (y1-y0)/(x1-x0);
    div_diff2 = (y2-y1)/(x2-x1);
    a = (div_diff2-div_diff1)/(x2-x0);
    
    b = div_diff2 + a*(x2-x1);
end

function [approx] = muller(x0, x1, x2, n)
    approx = [];
    copyx0 = x0;
    copyx1 = x1;
    copyx2 = x2;
    for i = 1:n
        x0 = copyx0;
        x1 = copyx1;
        x2 = copyx2;
        y0 = gx(x0, approx);
        y1 = gx(x1, approx);
        y2 = gx(x2, approx);
        while true
            [a, b, c] = get_coefficient(x0, x1, x2, y0, y1, y2);
            if b >= 0
                x3 = x2 - 2*c/(b+sqrt(b^2-4*a*c));
            else
                x3 = x2 - 2*c/(b-sqrt(b^2-4*a*c));
            end
            y3 = gx(x3, approx);
            if abs(x3-x2) <= 0.00005
                break;
            end
            x0 = x1;
            x1 = x2;
            x2 = x3;
            y0 = y1;
            y1 = y2;
            y2 = y3;
        end
        approx(i) = x2;
    end
end