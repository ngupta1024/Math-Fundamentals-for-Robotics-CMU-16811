clear;

data_file = '../cluttered_hallway.txt';
[t, error, avg_dist] = q4_e(data_file);

function [t, error, avg_dist] = q4_e(data_file)
    table_data = load(data_file);
    x = table_data(:, 1);
    y = table_data(:, 2);
    z = table_data(:, 3);
    
    t = [];
    avg_dist = [];
    for k = 4:-1:3
        [t_tmp, error, plane_data] = find_main_plane(x, y, z, k);
        t = [t, t_tmp];
        avg_dist = [avg_dist, get_avg_dist(x(plane_data), y(plane_data), z(plane_data), t_tmp)];
        x = x(plane_data==0);
        y = y(plane_data==0);
        z = z(plane_data==0);
        figure();
        scatter3(x, z, y, 'r.');
        set(gca,'Ydir','reverse');
        set(gca,'Zdir','reverse');
    end
    
    x2 = x(x<0);
    y2 = y(x<0);
    z2 = z(x<0);
    x1 = x(x>0);
    y1 = y(x>0);
    z1 = z(x>0);
    [~, tmp_count] = size(t);
    [t_tmp1, ~, ~] = find_main_plane(x1, y1, z1, 1.5);
    avg_dist(3) = get_avg_dist(x1, y1, z1, t_tmp1);
    t(:,3) = t_tmp1;
    while tmp_count < 4 || abs(t(1,3)/t(2,3)- t(1,4)/t(2,4)) > 0.5
        
        [t_tmp2, ~, ~] = find_main_plane(x2, y2, z2, 1.5);
        
        t(:,4) = t_tmp2;
        [~, tmp_count] = size(t);
    end
    avg_dist(4) = get_avg_dist(x2, y2, z2, t_tmp2);
    figure();
    scatter3(x2, z2, y2, 'r.');
    set(gca,'Ydir','reverse');
    set(gca,'Zdir','reverse');
    
    x = table_data(:, 1);
    y = table_data(:, 2);
    z = table_data(:, 3);
    plot_on_plane(x, y, z, t);
end

function [t, error, plane_data_ind] = find_main_plane(x, y, z, k)
    while true
        r = randi(length(x), 1, 3);
        r1 = r(1);
        r2 = r(2);
        r3 = r(3);
        xx = [x(r1); x(r2); x(r3)];
        yy = [y(r1); y(r2); y(r3)];
        zz = [z(r1); z(r2); z(r3)];
        [t, ~] = get_coefficient(xx, yy, zz);
    
        pred_z = t(1)*x+t(2)*y+t(3);
        dist = abs(pred_z-z)/sqrt(t(1)^2+t(2)^2+1);
        dist_th = dist <= 0.05;
        
        if sum(dist_th) > length(x)*(1/(k+1))
            break
        end
    end
    x = x(dist_th);
    y = y(dist_th);
    z = z(dist_th);
    plane_data_ind = dist <= 0.2;
    [t, error] = get_coefficient(x, y, z);
end

function [t, error] = get_coefficient(x, y, z)
    A = [x, y, ones(length(x), 1)];
    [u, s, vh] = svd(A);
    s_ = s';
    for k = 1:min(length(u), length(vh))
        if s_(k, k) ~= 0
            s_(k, k) = 1/s_(k, k);
        end
    end
    t = vh*s_*u'*z;
    
    zz = x*t(1) + y*t(2) + t(3);
    error = sqrt(sum((z - zz).^2));
end

function plot_on_plane(x, y, z, t)
    min_z = min(z);
    max_z = max(z);
    min_x = min(x);
    min_y = min(y);
    max_x = max(x);
    max_y = max(y);
    
    figure();
    scatter3(x, z, y, 'r.');
    set(gca,'Ydir','reverse');
    set(gca,'Zdir','reverse');

    figure();
    scatter3(x, z, y, 'r.');
    hold on;
    [~, plane_num] = size(t);
    for plane = 1:plane_num
        syms x_ y_ z_;
        planefunction=t(1, plane)*x_+t(2, plane)*y_+t(3, plane)-z_;
        zdep=solve(planefunction, y_);
        fsurf(zdep, 'FaceAlpha',0.5);
        hold on;
    end
    xlim([min_x-0.5 max_x+0.2]);
    zlim([min_y-0.2 max_y+0.2]);
    set(gca,'Ydir','reverse');
    set(gca,'Zdir','reverse');
end

function avg_dist = get_avg_dist(x, y, z, t)
    zz = abs(t(1)*x+t(2)*y+t(3)-z);
    dist = zz/sqrt(t(1)^2+t(2)^2+1);
    avg_dist = sum(dist)/length(zz);
end