
delta = 1;

DATA_train = csvread("training.csv",1);

phi = @(r) exp(-(r^2));


xs = DATA_train(:,1:end-1);
[nrow, ncol] = size(xs);
A = zeros(nrow,nrow);

for i=1:nrow
   for j=1:nrow
       r = norm(xs(i,:)-xs(j,:));
       A(i,j) = phi(r/delta);
   end
end



y = DATA_train(:,end);

c = A\y;


fnc =@(x) approximate(x, DATA_train(:,1:end-1), c, phi, delta);






DATA_val = csvread("validation.csv",1);
%RMSE = calcNormError(DATA_val, fnc);
%RMSE = calcNormError(DATA_train, fnc);


[n dim] = size(DATA_val);

err = 0;
for i=1:n
   fx = fnc(DATA_val(i,1:end-1));
   y = DATA_val(i,end);
   err = err + (fx-y)^2;
end

err = sqrt(err/n)
