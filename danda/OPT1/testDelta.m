function [ RMSE ] = testDelta( delta )
%TESTDELTA Summary of this function goes here
%   Detailed explanation goes here


DATA_train = csvread("training.csv",1);

phi = @(r) exp(-(r^2));

A = buildMatrix(DATA_train(:,1:end-1), phi, delta);
y = DATA_train(:,end);

c = A\y;


fnc =@(x) approximate(x, DATA_train(:,1:end-1), c, phi, delta);






DATA_val = csvread("validation.csv",1);
RMSE = calcNormError(DATA_val, fnc);
%RMSE = calcNormError(DATA_train, fnc);

end

