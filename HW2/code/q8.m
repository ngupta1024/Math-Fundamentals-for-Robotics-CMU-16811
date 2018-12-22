clear;

paths = load('../paths.txt');

%start_point = [0.8, 1.8];
%start_point = [2.2, 1.0];
start_point = [2.7, 1.4];
fire_center = [5, 5];

[interp_x, interp_y] = main(paths, start_point, fire_center);

function [up_path_x, up_path_y, low_path_x, low_path_y] = split_path(paths, fire_center)
    % This is used for path separation. One is the path group that above 
    % the fire ring, another path group include all paths below the fire ring.
    [m, n] = size(paths);
    x = paths(1:2:m, :);
    y = paths(2:2:m, :);
    up_y = find(y(:, fix(n/3)) > fire_center(2));
    low_y = find(y(:, fix(n/3)) < fire_center(2));
    up_path_x = x(up_y, :);
    up_path_y = y(up_y, :);
    low_path_x = x(low_y, :);
    low_path_y = y(low_y, :);
    
    figure(1);
    [m1, ~] = size(up_path_x);
    for i = 1:m1
        subplot(1,2,1);plot(up_path_x(i, :), up_path_y(i, :));
        hold on;
    end
    [m2, ~] = size(low_path_x);
    for j = 1:m2
        subplot(1,2,2);plot(low_path_x(j, :), low_path_y(j, :));
        hold on;
    end
    hold off;
end

function is_in = is_in_triangle(sp, p1, p2, p3)
    % This is to judge whether the start point is in the triangle (or on the line) 
    % formed by the given 3 points. If it is in the triangle 
    % or on the line formed by the 3 points, return true, 
    % else return false.
    is_in = true;
    A = [[p2(1)-p1(1), p3(1)-p1(1)];[p2(2)-p1(2), p3(2)-p1(2)]];
    b = [sp(1)-p1(1); sp(2)-p1(2)];
    if rank(A) < 2
        mins = min([p1; p2; p3], [], 1);
        maxs = max([p1; p2; p3], [], 1);
        if (sp(1)-mins(1))/(sp(2)-mins(2)) ~= (maxs(1)-mins(1))/(maxs(2)-mins(2))
            is_in = false;
        else
            if (sp(1)-mins(1))/(maxs(1)-mins(1)) < 0 || (sp(1)-mins(1))/(maxs(1)-mins(1)) > 1
                is_in = false;
            end
        end
        is_in = false;
    else
        v = A\b;
        if v(1) < 0 || v(1) > 1 || v(2) < 0 || v(2) > 1 || v(1)+v(2) > 1
            is_in = false;
        end
    end
end

function [path_set_x, path_set_y] = pick_path(up_path_x, up_path_y, low_path_x, low_path_y, start_point)
    % This is used to pick 3 paths from all the paths in order to 
    % construct the new path.
    up_center = [mean(up_path_x(:, 1)), mean(up_path_y(:, 1))];
    low_center = [mean(low_path_x(:, 1)), mean(low_path_y(:, 1))];
    if sum((start_point-up_center).^2) < sum((start_point-low_center).^2)
        path_x = up_path_x;
        path_y = up_path_y;
    else
        path_x = low_path_x;
        path_y = low_path_y;
    end
    
    [m, ~] = size(path_x);
    candidates = [];
    for i = 1:m
        for j = 2:m
            for k = 3:m
                p1 = [path_x(i, 1), path_y(i, 1)];
                p2 = [path_x(j, 1), path_y(j, 1)];
                p3 = [path_x(k, 1), path_y(k, 1)];
                is_in = is_in_triangle(start_point, p1, p2, p3);
                if is_in == true
                    candidates = [candidates; [i, j, k]];
                end
            end
        end
    end
    
    [~, n] = size(candidates);
    dist = [];
    for i = 1:n
        p1 = [path_x(candidates(i, 1), 1), path_y(candidates(i, 1), 1)];
        p2 = [path_x(candidates(i, 2), 1), path_y(candidates(i, 2), 1)];
        p3 = [path_x(candidates(i, 3), 1), path_y(candidates(i, 3), 1)];
        center = (p1+p2+p3)/3;
        dist = [dist; sum((center-start_point).^2)];
    end
    [~, ind] = min(dist);
    path_set_x = [path_x(candidates(ind, 1), :); 
                  path_x(candidates(ind, 2), :);
                  path_x(candidates(ind, 3), :)];
    path_set_y = [path_y(candidates(ind, 1), :); 
                  path_y(candidates(ind, 2), :);
                  path_y(candidates(ind, 3), :)];
end

function a = get_alpha(sp, p1, p2, p3)
    % This is to determine the weights alpha.
    P = [[p1(1)-p3(1), p2(1)-p3(1)];[p1(2)-p3(2), p2(2)-p3(2)]];
    b = [sp(1)-p3(1); sp(2)-p3(2)];
    a = P\b;
    a = [a; 1-a(1)-a(2)];
end

function [new_path_x, new_path_y] = get_new_path(a, path_set_x, path_set_y)
    % This is to generate the discrete points for new path
    new_path_x = a'*path_set_x;
    new_path_y = a'*path_set_y;
end

function interp_result = interpolate(data, t)
    interval_num = length(data)-1;
    data1 = data(fix(interval_num*t)+1);
    data2 = data(fix(interval_num*t)+2);
    interp_result = data1 + (data2-data1).*t;
end

function [interp_x, interp_y] = main(paths, start_point, fire_center)
    [up_path_x, up_path_y, low_path_x, low_path_y] = split_path(paths, fire_center);
    [path_set_x, path_set_y] = pick_path(up_path_x, up_path_y, low_path_x, low_path_y, start_point);
    
    p1 = [path_set_x(1, 1), path_set_y(1, 1)];
    p2 = [path_set_x(2, 1), path_set_y(2, 1)];
    p3 = [path_set_x(3, 1), path_set_y(3, 1)];
    a = get_alpha(start_point, p1, p2, p3);
    [new_path_x, new_path_y] = get_new_path(a, path_set_x, path_set_y);
    
    interp_x = interpolate(new_path_x, 0.01:0.01:0.99);
    interp_y = interpolate(new_path_y, 0.01:0.01:0.99); 
    
    figure(3);
    plot(start_point(1), start_point(2), 'b *', 'MarkerSize', 5);
    hold on;
    for i = 1:3
        plot(path_set_x(i, :), path_set_y(i, :), 'g');
        hold on;
    end
    plot(interp_x, interp_y, 'b', 'linewidth', 2);
    hold on;
    alpha=0:pi/90:2*pi;
    cir_x=1.5*cos(alpha)+5;
    cir_y=1.5*sin(alpha)+5;
    plot(cir_x, cir_y, 'm');
    legend('Starting point', 'selected path', 'selected path', 'selected path', 'interpolation result');
end
