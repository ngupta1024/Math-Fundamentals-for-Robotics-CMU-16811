

x = [0.7021, -0.0841];
%y = [0.0841, -0.7021];
y = [0,0];
xx = [-1:0.1:5];
yy = 0.*xx;

figure(1);
fimplicit(@(x,y) 2*x^2+2*y^2-1, [-1,5,-3,3]);
hold on;
fimplicit(@(x,y) x^2+y^2+2*x*y-x+y, [-1,5,-3,3]);
hold on;
plot(xx, yy, 'k :');
hold on;
plot(x,y,'r *', 'MarkerSize', 5);
hold off;
legend('2x^2+2y^2-1=0', 'x^2+y^2+2xy-x+y=0');
xlabel('x');
ylabel('y');