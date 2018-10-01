function hullEdges = listEdges(hullFaces,X)

% @pre    : X(i,j) is the j^th coordinate of the node with index i
% @pre    : hullFaces(i,j) is the j^th node of the face with index i on the convex hull of X
% @return : hullEdges(i,:) describes the edge with index i ; its six components are the indexes of : 
% - the node at its origin
% - the node at its tip
% - the face on its left
% - the face on its right
% - the coordinates of the normalized vector aligned with the edge
% a face with index 0 means there is no face on that side
%
% MELCHIOR Samuel - UCL, 2008-09.

nEdges = 0;
hullEdges = zeros(0,7);
for k=1:3*size(hullFaces,1)
    if mod(k,3) == 1,
        face = hullFaces(ceil(k/3),:);
        normal = -normalToFaces(face,X);
    end
    newEdge = sort([face(1+mod(k-1,3)) face(1+mod(k,3))]);
    [idx,insert] = findRow(hullEdges(1:nEdges,1:2),newEdge);
    if idx > 0,
        isRight = hullEdges(idx,4) == 0;
        hullEdges(idx,3 + isRight) = ceil(k/3);
    else,
        if insert > 0,
            vecEdge = X(newEdge(2),:)-X(newEdge(1),:);
            hullEdges(insert:nEdges+1,:) = [newEdge 0 0 vecEdge/norm(vecEdge);hullEdges(insert:nEdges,:)];
            isRight = (crossNorm(normal,vecEdge)*(X(face(1+mod(k+1,3)),:) - X(face(1+mod(k-1,3)),:))') < 0;
            hullEdges(insert,3 + isRight) = ceil(k/3);
            nEdges = nEdges+1;
        else,
            assert(false);
        end
    end
end

function [idx, insert] = findRow(M, row)

if (isempty(M)),
    idx = 0;
    insert = 1;
else,
    m = size(M, 2);
    low = 1;
    high = size(M, 1) + 1;
    while (low ~= high),
        mid = floor((low + high) / 2);
        for i=1:m+1,
            if (i == m + 1),
                idx = mid;
                insert = 0;
                return;
            end
            if row(i) < M(mid, i),
                high = mid;
                break;
            elseif row(i) > M(mid, i),
                low = mid + 1;
                break;
            end
        end
    end

    idx = 0;
    insert = low;
end

