function [B, index_kept] = prune_library2(A,min_angle)
%        B = prune_library(A,min_angle)
%
% remove columns from A such that the minimum angle between the columns
% of B in no smaller than max_angle
%
% Author: Jose Bioucas Dias. June, 2011
%

[L,m] = size(A);  % L = number of bands; m = number of materilas
% na=sqrt(sum((Y_5)'.^2));
% nb=sqrt(sum((Y_6)'.^2));
% A_norm=Y_5'./repmat(na,224,1);
% B_norm=Y_6'./repmat(nb,224,1);
%normalize A
nA = sqrt(sum(A.^2));
A_norm = A./repmat(nA,L,1);
% compute angles
angles = abs(acos(A_norm'*A_norm))*180/pi;
% discard vectors with angles less than min_angle
index = 1;
index_kept = [];

for i=1:m
   if angles(i,i) ~= inf
       B(:,index) = A(:,i);
       index_kept = [index_kept i];
       angles(:,angles(i,:) < min_angle ) = inf;
       index = index + 1;
   end

end





