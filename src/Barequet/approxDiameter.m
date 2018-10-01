function st = approxDiameter(X)

% @pre    : X(i,j) is the j^th coordinate of the node with index i [j=1:d]
% @return : st is a 1/sqrt(d)-approximation of the diameter of X
%
% MELCHIOR Samuel - UCL, 2008-09.

[X_Min,nodeMin] = min(X);
[X_Max,nodeMax] = max(X);
[l,ind_d] = max(X_Max-X_Min);
st = X(nodeMax(ind_d),:) - X(nodeMin(ind_d),:);
