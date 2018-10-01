function c = crossNorm(a,b)

% @pre    : a and b \in |R^3
% @return : c is a x b
%
% MELCHIOR Samuel - UCL, 2008-09.

c = [a(2)*b(3)-a(3)*b(2) a(3)*b(1)-a(1)*b(3) a(1)*b(2)-a(2)*b(1)];
c = sqrt(c*c')\c; 