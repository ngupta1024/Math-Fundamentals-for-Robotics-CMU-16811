clear;

plot_equation();
[x_low] = zero_solution(10.9);
[x_high] = zero_solution(14);

function plot_equation()
    interval1 = [10 11];
    figure(1);
    fplot(@(x) tan(x)-x, interval1);
    hold on;
    fplot(@(x) 0, interval1)
    hold off;
    
    interval2 = [13 15];
    figure(2);
    fplot(@(x) tan(x)-x, interval2);
    hold on;
    fplot(@(x) 0, interval2)
    hold off;
end

function [y] = fx(x)
    y = tan(x) - x;
end

function [y] = derivative_fx(x)
    y = (1/cos(x))^2 - 1;
end

function [x_new] = newton_method(x)
    x_new = x - fx(x)/derivative_fx(x);
end

function [x] = zero_solution(x0)
    x = x0;
    while abs(fx(x)) > 0.000000001
        x = newton_method(x);
    end
end