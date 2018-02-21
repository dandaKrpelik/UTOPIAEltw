function [ err ] = calcNormError( data, fnc )
%CALCNORMERROR Summary of this function goes here
%   Detailed explanation goes here

[n dim] = size(data);

err = 0;
for i=1:n
   fx = fnc(data(i,1:end-1));
   y = data(i,end);
   
   err = err + (fx-y)^2;
end

err = sqrt(err/n);

end

