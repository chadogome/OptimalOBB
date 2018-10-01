function theta1 = theta2toTheta1(theta2,phi)

c2 = cos(theta2).^2;
s2 = sin(phi).^2;
% c2 = s2 (1 - x^2)/(s2 + (1-s2)*x^2) // x = cos(theta1)
% c2*(s2 + (1-s2)*x^2) = s2 (1 - x^2)
% (c2+s2*(1-c2))x^2 = s2*(1-c2)
theta1 = -real(acos(-sign(theta2).*sqrt(s2*(1 - c2)./(c2*(1-s2) + s2))));
[sin(theta1) -cos(theta2).*sqrt(cos(theta1).^2/s2+sin(theta1).^2)]