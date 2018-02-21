clear all

n = 50;
d = linspace(0.001,5,n);
y = zeros(n,1);

for i=1:n
    y(i) = testDelta(d(i));
end

plot(d,y)