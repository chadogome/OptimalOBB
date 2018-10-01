function [V, R] = rotatingCalipersAround(data, axis)

if (axis(1) == 0 && axis(2) == 0),
    Rt = eye(3);
else
    xaxis = axis(2)/sqrt(axis(1)*axis(1)+axis(2)*axis(2));
    yaxis = -axis(1)/sqrt(axis(1)*axis(1)+axis(2)*axis(2));
    cosangle = axis(3)/sqrt(axis(1)*axis(1)+axis(2)*axis(2)+axis(3)*axis(3));
    sinangle = sqrt(1-cosangle*cosangle);
    Rt = [xaxis*xaxis*(1-cosangle) + cosangle,  xaxis*yaxis*(1-cosangle),            -yaxis*sinangle ;
        xaxis*yaxis*(1-cosangle),             yaxis*yaxis*(1-cosangle) + cosangle, xaxis*sinangle  ;
        yaxis*sinangle,                       -xaxis*sinangle,                     cosangle       ];
end
rotatedData = data*Rt;
[vol, angle] = rotatingCalipers(rotatedData(:, 1:2));
V = vol*(max(rotatedData(:, 3)) - min(rotatedData(:, 3)));
R = [cos(angle) sin(angle) 0 ; -sin(angle) cos(angle) 0 ; 0 0 1]*Rt';