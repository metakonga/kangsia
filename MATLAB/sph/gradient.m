function gradient(x, y1, y2, div1, div2)
len = length(x);
grad = zeros(2,6);
for i = 1:len
    grad(1,i) = (y1(i+1) - y1(i)) / div1;
    grad(2,i) = (y2(i+1) - y2(i)) / div2;
end

plot(x, grad(1,:));
hold on
plot(x, grad(2,:));