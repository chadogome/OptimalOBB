function R_n = computeR_n(theta1,phi,R_x,normal2)

% @pre    : theta1 is the angle between R_n(1,:) and R_x(1,:) around R_x(3,:)
% @pre    : phi is the angle R_x(2,:) and R_x(3,:) around R_x(1,:)
% @pre    : R_x(i,:) is the frame determined by the two edges of the current iteration
% @return : R_n(i,:) is the i^th axis for the frame aligned with the bounding box
%
% MELCHIOR Samuel - UCL, 2008-09.
% 

ct = cos(theta1); st = sin(theta1);
cp = cos(phi); sp = sin(phi);
n1_R_x = [ct, st, 0]*R_x;
n2_R_x = ([-st*sp, ct*sp, -ct*cp]*R_x)/sqrt(ct*ct + st*st*sp*sp);
if n2_R_x*normal2' < 0,
    n2_R_x = -n2_R_x;
end
n3_R_x = [n1_R_x(2)*n2_R_x(3)-n1_R_x(3)*n2_R_x(2) n1_R_x(3)*n2_R_x(1)-n1_R_x(1)*n2_R_x(3) n1_R_x(1)*n2_R_x(2)-n1_R_x(2)*n2_R_x(1)];
R_n = [n1_R_x;n2_R_x;n3_R_x];
