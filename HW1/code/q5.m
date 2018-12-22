clear;

P=[[0;0;-1],[-1;-1;-3],[-1;-1;0],[-2;-3;-1],[-2;-1;-1],[-0.5;-2;-1.5],[-0.2;-0.5;-1],[-1.3;-0.3;0],[-1.5;-1.5;-1.5],[-0.4;-0.5;0.8]];
Q=[[2.5;0;2],[3.5;-1;4],[3.5;1.5;2.5],[4.5;1.6;4.6],[4.5;0.6;2.8],[3;0.7;4],[2.7;0.4;2.4],[3.8;0.8;1.8],[4;0.5;3.5],[3;2;1.5]];

[A,t]=transform_info(P,Q);
plot_points(P,Q,A,t);

function [A,t]=transform_info(P,Q)
% Inputs:
%   P: The input 3*n matrix where each column is the 3D coordinates
%      of a point.
%   Q: The input 3*n matrix where each column is the 3D coordinates 
%      anslated.
% Outputs:
%   A: 3*3 matrix that performs the rotation
%   t: 3*1 vector that performs the translation
    [~,n]=size(P);
    % get the centers of p and q
    p_avg=sum(P,2)/n;
    q_avg=sum(Q,2)/n;

    % To calculate sum((q-q_avg)*((p-p_avg)^T))/n
    B=zeros(3,3);
    for k=1:n
        B=B+(Q(:,k)-q_avg)*((P(:,k)-p_avg)');
    end
    B=B/n;
    % perform the SVD decomposition to B
    [U,S,V]=svd(B);
    % calculate A and t
    A=U*V';
    t=q_avg-A*p_avg;
end

function plot_points(P,Q,A,t)
% Inputs:
%   P: The input 3*n matrix where each column is the 3D coordinates
%      of a point.
%   Q: The input 3*n matrix where each column is the 3D coordinates 
%      of the corresponding point in P after translated.
%   A: 3*3 matrix that performs the rotation
%   b: 3*1 vector that performs the translation
% Outputs:  
    [~,n]=size(P);
    % get the centers of p and q
    p_avg=sum(P,2)/n;
    q_avg=sum(Q,2)/n;
    % Rotate and translate the points in P using the obtained 
    % rotation and transition matrix A and t
    R=A*P+t;
    % Rotate and translate the center point of P
    r_center=A*p_avg+t;

    % plot the original points and the center points
    plot3(P(1,:),P(2,:),P(3,:),'b.','MarkerSize',15);
    hold on;
    plot3(p_avg(1),p_avg(2),p_avg(3),'b.','Marker','*','MarkerSize',5);
    hold on;
    plot3(Q(1,:),Q(2,:),Q(3,:),'r.','MarkerSize',15);
    hold on;
    plot3(q_avg(1),q_avg(2),q_avg(3),'r.','Marker','*','MarkerSize',5);
    hold on;
    % plot the points after rotation and translation by the obtained A and t
    plot3(R(1,:),R(2,:),R(3,:),'g.','MarkerSize',15);
    hold on;
    plot3(r_center(1),r_center(2),r_center(3),'g.','Marker','o','MarkerSize',10);
    grid on;
    legend('points P', 'center point of P', 'points Q', 'center point of Q', 'points R=AP+t', 'center point of R');
end

