function hullNodes = nodeToFaces(hullFaces,X)

% @pre    : X(i,j) is the j^th coordinate of the node with index i
% @pre    : hullFaces(i,j) is the j^th node of the face with index i on the convex hull of X
% @return : hullNodes(i,:) is a list of the faces incident to the node with index i, sorted counterclockwise
%
% MELCHIOR Samuel - UCL, 2008-09.

hullNodes = cell( max(hullFaces(:)),1);
for k=1:3*size(hullFaces,1)
    if mod(k,3) == 1
        face = hullFaces(ceil(k/3),:);
        center = mean(X(face,:));
        normal = -normalToFaces(face,X);
    end
    newNode = face(1+mod(k-1,3));
    faces = hullNodes{newNode};
    l = [];
    for l=1:length(faces)
        center_l = mean(X(hullFaces(faces(l),:),:));
        if crossNorm(center-X(newNode,:),center_l-X(newNode,:))*mean([normal;-normalToFaces(hullFaces(faces(l),:),X)])'> 0 %normal' < 0
            l = l-1;
            break;
        end
    end
    if isempty(l)
        l = 1;
    else
        if l==0
            for l=length(faces):-1:1
                center_l = mean(X(hullFaces(faces(l),:),:));
                if crossNorm(center-X(newNode,:),center_l-X(newNode,:))*mean([normal;-normalToFaces(hullFaces(faces(l),:),X)])'< 0 %normal' < 0
                    l = l+1;
                    break;
                end
            end
            l = l-1;
        end
        l = l+1;
    end
    hullNodes{newNode}(l:end+1) = [ceil(k/3) faces(l:end)];
end
