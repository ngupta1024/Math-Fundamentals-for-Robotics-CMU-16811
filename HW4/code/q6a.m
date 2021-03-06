%% create obstacle field
close all
clear all
waypoints=300;
N=101;
OBST = [20,30;60,40;70,85];
epsilon = [25; 20; 30];

obs_cost = double(zeros(N));
for i=1:size(OBST,1)

    t = zeros(N);
    t(OBST(i,1),OBST(i,2)) = 1; %point obstacles
    
    t_cost = double(bwdist(t));
    t_cost(t_cost>epsilon(i))=epsilon(i);
    t_cost = 1/(2*epsilon(i))*(t_cost-epsilon(i)).^2;
    
    obs_cost = obs_cost + t_cost(1:N, 1:N);
end


% obstacle cost gradient
gx = diff(double(obs_cost),1,1);
gy = diff(double(obs_cost),1,2);


figure(1);
surf(1:N,1:N,double(obs_cost')), view(-25,70);
xlabel('X')
xlim([0 100]);
ylim([0 100]);
hold on;

%% initial path
%world params
SX = 10; % START
SY = 10;
GX = 90; % GOAL
GY = 90;


traj = zeros(2,waypoints);
traj(1,1)=SX;
traj(2,1)=SY;
dist_x = GX-SX;
dist_y = GY-SY;
for i=2:waypoints
    traj(1,i)=traj(1,i-1)+dist_x/(waypoints-1);
    traj(2,i)=traj(2,i-1)+dist_y/(waypoints-1);
end

path_init = traj';
tt=size(path_init,1);
path_init_values = zeros(size(path_init,1),1);
for i=1:tt
    path_init_values(i)=obs_cost(floor(path_init(i,1)),floor(path_init(i,2)));
end
plot3(path_init(:,1),path_init(:,2),path_init_values,'.r','MarkerSize',15);
hold on

path = path_init;


%% Optimize it...
% your code comes here
for i=2:tt-1
    x = path(i, 1);
    y = path(i, 2);
    interp_x1 = (x - fix(x))*(gx(fix(x)+1, fix(y)) - gx(fix(x), fix(y)))+ gx(fix(x), fix(y));
    interp_x2 = (x - fix(x))*(gx(fix(x)+1, fix(y)+1) - gx(fix(x), fix(y)+1))+ gx(fix(x), fix(y)+1);
    grad_x = (y - fix(y))*(interp_x2-interp_x1) + interp_x1;
    interp_y1 = (x - fix(x))*(gy(fix(x)+1, fix(y)) - gy(fix(x), fix(y)))+ gy(fix(x), fix(y));
    interp_y2 = (x - fix(x))*(gy(fix(x)+1, fix(y)+1) - gy(fix(x), fix(y)+1))+ gy(fix(x), fix(y)+1);
    grad_y = (y - fix(y))*(interp_y2-interp_y1) + interp_y1;
    
    linex = gx(:, fix(y));
    liney = gy(fix(x), :);
    if grad_x > 0
        stepx = sum(linex>0);
    else
        stepx = sum(linex<0);
    end
    if grad_y > 0
        stepy = sum(liney>0);
    else
        stepy = sum(liney<0);
    end
    
    scale = 0.1;
    path(i, 1) = x - scale*grad_x;
    path(i, 2) = y - scale*grad_y;
end

figure(2);
surf(1:N,1:N,double(obs_cost')), view(-25,70);
xlabel('X')
xlim([0 100]);
ylim([0 100]);
hold on;
path_values = zeros(tt,1);
for i=1:tt
    path_values(i)=obs_cost(floor(path(i,1)),floor(path(i,2)));
end
figure(2)
hold on;
plot3(path_init(:,1),path_init(:,2),path_init_values,'.r','MarkerSize',15);
hold on;
plot3(path(:,1),path(:,2),path_values,'.g','MarkerSize',25);
title('One iteration')
hold off;

for i=2:tt-1
    while true
        x = path(i, 1);
        y = path(i, 2);
        interp_x1 = (x - fix(x))*(gx(fix(x)+1, fix(y)) - gx(fix(x), fix(y)))+ gx(fix(x), fix(y));
        interp_x2 = (x - fix(x))*(gx(fix(x)+1, fix(y)+1) - gx(fix(x), fix(y)+1))+ gx(fix(x), fix(y)+1);
        grad_x = (y - fix(y))*(interp_x2-interp_x1) + interp_x1;
        interp_y1 = (x - fix(x))*(gy(fix(x)+1, fix(y)) - gy(fix(x), fix(y)))+ gy(fix(x), fix(y));
        interp_y2 = (x - fix(x))*(gy(fix(x)+1, fix(y)+1) - gy(fix(x), fix(y)+1))+ gy(fix(x), fix(y)+1);
        grad_y = (y - fix(y))*(interp_y2-interp_y1) + interp_y1;
        if (grad_x^2+grad_y^2) < 0.00001
            break
        end
    
        linex = gx(:, fix(y));
        liney = gy(fix(x), :);
        if grad_x > 0
            stepx = sum(linex>0);
        else
            stepx = sum(linex<0);
        end
        if grad_y > 0
            stepy = sum(liney>0);
        else
            stepy = sum(liney<0);
        end
    
        scale = 0.1;
        path(i, 1) = x - scale*grad_x;
        path(i, 2) = y - scale*grad_y;
    end
end

%% plot the trajectories
path_values = zeros(tt,1);
for i=1:tt
    path_values(i)=obs_cost(floor(path(i,1)),floor(path(i,2)));
end
figure(1)
hold on;
plot3(path(:,1),path(:,2),path_values,'.g','MarkerSize',25);

hold off;
