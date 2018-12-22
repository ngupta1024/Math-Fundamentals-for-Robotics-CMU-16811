clear;

x1 = 2.25;
data_x1 = [1.0,1.5,2.0,2.5,3.0,3.5,4.0];
data_y1 = get_fx_1(data_x1);
coefficient1 = interpolate(data_x1, data_y1);
y1 = estimate_y(coefficient1, x1, data_x1);

x2 = 0.05;
n = 40;
data_x2 = get_x(n);
data_y2 = get_fx_2(data_x2);
coefficient2 = interpolate(data_x2, data_y2);
y2 = estimate_y(coefficient2, x2, data_x2);

n_list = [2,4,6,8,10,12,14,16,18,20,40];
En = error_estimate(n_list);

function [data_y] = get_fx_1(data_x)
    [~, n] = size(data_x);
    data_y = zeros(1, n);
    for i = 1:n
        data_y(i) = (log(data_x(i))/log(6))^1.5;
    end
end

function [data_x] = get_x(n)
    data_x = zeros(1, n+1);
    for i = 0:n
        data_x(i+1) = i*2/n-1;
    end
end

function [data_y] = get_fx_2(data_x)
    [~, n] = size(data_x);
    data_y = zeros(1, n);
    for i = 1:n
        data_y(i) = 6/(1+25*(data_x(i)^2));
    end
end

function coefficient = interpolate(data_x, data_y)
    [~, n] = size(data_x);
    div_diff_table = zeros(n, n);
    div_diff = data_y;
    div_diff_table(:, 1) = div_diff;
    for i = 2:n
        next_div_diff = zeros(1, n);
        for j = i:n
            next_div_diff(j) = (div_diff(j)-div_diff(j-1))/(data_x(j)-data_x(j-i+1));
        end
        div_diff = next_div_diff;
        div_diff_table(:, i) = div_diff;
    end
    coefficient = diag(div_diff_table);
end

function y = estimate_y(coefficient, x, data_x)
    y = 0;
    for i = 1:length(coefficient)
        item = 1;
        for j = 1:i-1
            item = item*(x-data_x(j));
        end
        y = y + coefficient(i)*item;
    end
end

function [En] = error_estimate(n_list)
    discretize = 1000;
    x = -1:(2/discretize):1;
    y = get_fx_2(x);
    [~, l] = size(n_list);
    En = zeros(1, l);
    for i = 1:l
        data_x = get_x(n_list(i));
        data_y = get_fx_2(data_x);
        esti_y = zeros(1, discretize+1);
        for j = 1:discretize+1
            coeff = interpolate(data_x, data_y);
            esti_y(j) = estimate_y(coeff, x(j), data_x);
        end
        En(i) = max(abs(esti_y - y));
    end
end