clear;

[avg_error, avg_error2, avg_error3] = main();

function [avg_error, avg_error2, avg_error3] = main()
    x0 = 2;
    y0 = 1;
    h = -0.05;
    x = 2:h:1;
    true_y = get_true_value(x);
    
    euler_y = euler(x, x0, y0, h);
    avg_error = sum(abs(euler_y - true_y)) / length(x);
    
    rk4_y = RK4(x, x0, y0, h);
    avg_error2 = sum(abs(rk4_y - true_y)) / length(x);
    
    ab4_y = AB4(x, x0, y0, h);
    avg_error3 = sum(abs(ab4_y - true_y)) / length(x);
    
    plot_values(x, true_y, euler_y);
    plot_values(x, true_y, rk4_y);
    plot_values(x, true_y, ab4_y);
end

function true_y = get_true_value(x)
    true_y = (x-1).^(1/3);
end

function fn = f(x)
    fn = 1/(3*x^2);
end

function estimate_y = euler(x, x0, y0, h)
    estimate_y = [y0];
    pre_y = y0;
    for i = 2:length(x)
        y_new = pre_y + h*(f(pre_y));
        estimate_y = [estimate_y, y_new];
        pre_y = y_new;
    end
end

function estimate_y = RK4(x, x0, y0, h)
    estimate_y = [y0];
    pre_y = y0;
    for i = 2:length(x)
        k1 = h*f(pre_y);
        k2 = h*f(pre_y+k1/2);
        k3 = h*f(pre_y+k2/2);
        k4 = h*f(pre_y+k3);
        y_new = pre_y + 1/6*(k1+2*k2+2*k3+k4);
        estimate_y = [estimate_y, y_new];
        pre_y = y_new;
    end
end

function estimate_y = AB4(x, x0, y0, h)
    estimate_y = [y0];
    pre_y = y0;
    f4 = f(1.04768955317165);
    f3 = f(1.03228011545637);
    f2 = f(1.01639635681485);
    f1 = f(1);
    for i = 2:length(x)
        y_new = pre_y + h/24*(55*f1-59*f2+37*f3-9*f4);
        estimate_y = [estimate_y, y_new];
        f4 = f3;
        f3 = f2;
        f2 = f1;
        f1 = f(y_new);
        pre_y = y_new;
    end
end

function plot_values(x, true_y, estimate_y)
    fig = figure();
    data = [x; true_y; estimate_y; true_y-estimate_y; abs(true_y-estimate_y)]';
    colnames = {'Xi', 'True value', 'Estimate value', 'True-Estimate', 'abs(True-Estimate)'};
    t = uitable(fig, data, colnames, 'Position',[0 0 500 500]);
    figure();
    plot(x, true_y, 'r');
    hold on;
    plot(x, estimate_y, 'g');
    legend('True value', 'Estimate value');
end