function [Ropt, Vopt] = pca(points, params)

data_means = mean(points);
B = points - ones(size(points, 1), 1)*data_means;
[eigvec, eigval] = eig(B'*B);

if strcmp(params.heuristic,'min')
    [foo, index] = min(diag(eigval));
elseif strcmp(params.heuristic,'max')
    [foo, index] = max(diag(eigval));
else
    Ropt = eigvec';
    rotated = points*eigvec;
    Vopt = prod(max(rotated) - min(rotated));
end

if strcmp(params.heuristic,'min') || strcmp(params.heuristic,'max')
    [Vopt, Ropt] = rotatingCalipersAround(points, eigvec(:, index(1)));
    
%     if (axis(1) == 0 && axis(2) == 0),
%         Rt = eye(3);
%     else
%         xaxis = axis(2)/sqrt(axis(1)*axis(1)+axis(2)*axis(2));
%         yaxis = -axis(1)/sqrt(axis(1)*axis(1)+axis(2)*axis(2));
%         cosangle = axis(3)/sqrt(axis(1)*axis(1)+axis(2)*axis(2)+axis(3)*axis(3));
%         sinangle = sqrt(1-cosangle*cosangle);
%         Rt = [xaxis*xaxis*(1-cosangle) + cosangle,  xaxis*yaxis*(1-cosangle),            -yaxis*sinangle ;
%             xaxis*yaxis*(1-cosangle),             yaxis*yaxis*(1-cosangle) + cosangle, xaxis*sinangle  ;
%             yaxis*sinangle,                       -xaxis*sinangle,                     cosangle       ];
%     end
%     rotatedX = points*Rt;
%     [vol, angle] = rotatingCalipers(rotatedX(:, 1:2));
%     Vopt = vol*(max(rotatedX(:, 3)) - min(rotatedX(:, 3)));
%     Ropt = [cos(angle) sin(angle) 0 ; -sin(angle) cos(angle) 0 ; 0 0 1]*Rt';
end
