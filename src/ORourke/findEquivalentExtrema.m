function nodes = findEquivalentExtrema(normal,nodes,allFaces,allNodes,X,tol)

% @pre    : 
% @pre    : hull***** is defined below [optional]
% @pre    : only the pair of edges i0 and j0(<i0) is tested [optional]
% @return : R_opt(i,:) is the i^th axis for the frame aligned with the minimum volume oriented bounding box
% @return : V_opt is the minimal volume over all the oriented bounding box
%
% MELCHIOR Samuel - UCL, 2008-09.

faces = allNodes{nodes(end)};
for i=1:length(faces);
    for j=1:3
        newNode = allFaces(faces(i),j);
        if isempty(find(nodes == newNode,1)) && (X(newNode,:)*normal' >= X(nodes(1),:)*normal'-tol*max(abs(X(nodes(1),:)*normal'),1))
            nodes = findEquivalentExtrema(normal,[nodes newNode],allFaces,allNodes,X,tol);
        end
    end
end