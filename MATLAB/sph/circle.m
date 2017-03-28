function circle(a,b,r)

t = linspace(0, 2*pi, 1000); % À§»ó °ª

y_r=r*cos(t);

y_i=r*sin(t);




y_r2=y_r+a;

y_i2=y_i+b;

figure, plot(y_r2,y_i2) 
