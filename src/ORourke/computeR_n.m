function R_n = computeR_n(theta1,phi,R_x,normal2)

% @pre    : theta1 is the angle between R_n(1,:) and R_x(1,:) around R_x(3,:)
% @pre    : phi is the angle R_x(2,:) and R_x(3,:) around R_x(1,:)
% @pre    : R_x(i,:) is the frame determined by the two edges of the current iteration
% @return : R_n(i,:) is the i^th axis for the frame aligned with the bounding box
%
% MELCHIOR Samuel - UCL, 2008-09.
% 

n1 = [cos(theta1) sin(theta1) 0];
n2 = [-sin(theta1)*sin(phi) cos(theta1)*sin(phi) -cos(theta1)*cos(phi)]/sqrt(cos(theta1)^2+sin(theta1)^2*sin(phi)^2);
sigma = sign(n2*R_x*normal2');
if sigma == 0
    sigma = 1;
end
n2 = n2*sigma;
%n3 = cross(n1,n2); %[n2(3)*sin(theta1) -n2(3)*cos(theta1) -n2(1)/sin(theta1)];   % !!! it bugged for theta1 = 0
n3 = [n1(2)*n2(3) -n1(1)*n2(3) n1(1)*n2(2)-n1(2)*n2(1)];
R_n = [n1;n2;n3]*R_x;