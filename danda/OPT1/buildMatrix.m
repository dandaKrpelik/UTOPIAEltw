function [ A ] = buildMatrix( xs , phi, delta )
%BUILDMATRIX Summary of this function goes here
%   Detailed explanation goes here

[nrow, ncol] = size(xs);

A = zeros(nrow,nrow);

for i=1:nrow
   for j=1:nrow
       r = norm(xs(i,:)-xs(j,:));
       
       A(i,j) = phi(r/delta);
   end
end



end

