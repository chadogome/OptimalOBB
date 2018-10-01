function c = normalToFaces(faces,X)

% @pre    : faces(i,j) is the j^th node of the i^th face
% @pre    : X(i,j) is the j^th coordinate of the node with index i
% @return : c(i,j) is the j^th compotent of the normal to the i^th face
%
% MELCHIOR Samuel - UCL, 2008-09.

a = X(faces(:,2),:)-X(faces(:,1),:);
b = X(faces(:,3),:)-X(faces(:,1),:);

c = [a(:,2).*b(:,3)-a(:,3).*b(:,2) a(:,3).*b(:,1)-a(:,1).*b(:,3) a(:,1).*b(:,2)-a(:,2).*b(:,1)];
c = sqrt(diag(diag(c*c')))\c; 