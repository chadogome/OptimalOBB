function [Ropt,Vopt] = barequetFast(X,d)

%[Ropt,Vopt] = barequetFast(X,d)
% @pre    : X(i,j) is the j^th coordinate of the node with index i
% @pre    : d (= c/epsilon) determines the number of axis around which we search the minimum volume oriented bounding box with rotating calipers
% @post   : the oriented bounding box is an (1+epsilon)-approximation of the minimum volume oriented bounding box (assuming c is well chosen)
% @return : R_opt(i,:) is the i^th axis for the frame aligned with the oriented bounding box
% @return : V_opt is the  volume of the oriented bounding box


% Computation of a 15-approximation of the bounding box
st = approxDiameter(X);
% R = zeros(3);
% R(1,:) = st/norm(st);
% Q = X - (X*R(1,:)')*R(1,:);
% s_t_ = approxDiameter(Q);
% R(2,:) = s_t_/norm(s_t_);
% R(3,:) = cross(R(1,:),R(2,:));

[v, R] = rotatingCalipersAround(X, st);
S = X*R';
bbox = max(S)-min(S);

% Loop on all the directions induced by the grid define by bbox and c_epsilon
Vopt = Inf;
for i=-d:d
    for j=-d:d
        for k=0:d
            v = bbox.*[i j k];
            [V, RR] = rotatingCalipersAround(X, v);
            if V < Vopt
                Vopt = V;
                Ropt = RR;
            end
        end
    end
end