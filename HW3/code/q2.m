clear;

x = 0:0.01:1;
y = load('../problem2.txt');

figure(1);
plot(x, y);

y0 = y(1);
y1 = y(length(y));

g = (y1-y0)*x + y0;
figure(2);
plot(x, y);
hold on;
plot(x, g);

t = y-g;
figure(3);
plot(x, t);

base1 = x;
base2 = sin(6*pi*x);
[A, error] = get_coefficient(y, x, [base1;base2]);
yy = x*A(1) + A(2)*sin(6*pi*x);
figure(4);
plot(x, yy, 'linewidth', 2);

function [A, error] = get_coefficient(y, x, base)
    base1 = base(1, :);
    base2 = base(2, :);
    
    M = [base1; base2]';
    [u, s, vh] = svd(M);
    s_ = zeros(length(vh), length(u));
    for m = 1:length(vh)
        if s(m, m) ~= 0
            s_(m, m) = 1/s(m, m);
        end
    end
    A = vh*s_*u'*y';
    yy = A(1)*base1 + A(2)*base2;
    error = sqrt(sum((y - yy).^2));
end