k1 = 0.5;
k2 = 2;

x = -2:0.1:2;
y1 = 0.5*k1*x.^2;
y2 = 0.5*k2*x.^2;
y3 = gauss(x,1/sqrt(k1));
y4 = gauss(x,1/sqrt(k2));

plot(x,y1,x,y2,x,y3,x,y4)
legend('0.5','2','0.5','2')

function y = gauss(x,s2)
    y = 1/sqrt(2*pi*s2)*exp(-(x).^2/(2*s2));
end 