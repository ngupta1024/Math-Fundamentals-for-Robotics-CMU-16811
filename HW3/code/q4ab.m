clear;

data_file1 = '../clear_table.txt';
[t1, error1, evg_dist1] = q4_a(data_file1);

data_file2 = '../cluttered_table.txt';
[t2, error2, evg_dist2] = q4_a(data_file2);

function [t, error, avg_dist] = q4_a(data_file)
    table_data = load(data_file);
    x = table_data(:, 1);
    y = table_data(:, 2);
    z = table_data(:, 3);

    [t, error] = get_coefficient(x, y, z);
    plot_on_plane(x, y, z, t);
    avg_dist = get_avg_dist(x, y, z, t);
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
    syms x_ y_ z_;
    planefunction=t(1)*x_+t(2)*y_+t(3)-z_;
    zdep=solve(planefunction, y_);
    fsurf(zdep,[min_x-0.1,max_x+0.1,min_z-0.1,max_z+0.1], 'FaceAlpha',0.5);
    set(gca,'Ydir','reverse');
    set(gca,'Zdir','reverse');
    zlim([min_y-0.1 max_y+0.1]);
end

function avg_dist = get_avg_dist(x, y, z, t)
    zz = abs(t(1)*x+t(2)*y+t(3)-z);
    dist = zz/sqrt(t(1)^2+t(2)^2+1);
    avg_dist = sum(dist)/length(zz);
end