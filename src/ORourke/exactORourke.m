function [Vopt,Ropt] = exactORourke(Xname)

if ((nargin == 1) && ischar(Xname) && strcmp(Xname, 'VersionCheck42')),
    disp('O''Rourke v7m SS Matlab');
    disp('C.-T. Chang, B. Gorissen, S. Melchior*');
    Vopt = [];
    Ropt = [];
    return;
end

if ((nargin < 1) || ~ischar(Xname) || (exist(Xname, 'file') ~= 2)),
    disp('Argument of this Matlab function must be a string corresponding to a .mat file (e.g. argument ''input.mat'').');
    disp('This .mat file should contain a variable X of size Nx3 describing the coordinates of the set of N points.');
    disp(' ');
    disp('Return values Vopt and Ropt will be stored in the file ''output.mat''.');
    disp('Vopt = Minimal OBB volume of the set of points (tolerance = 1e-9).');
    disp('Ropt = Rotation matrix associated to the optimal OBB.');
    disp(' ');
    disp('Bugs and errors can be reported on https://github.com/chagome/OptimalOBB/issues.');
    Vopt = [];
    Ropt = [];
    return;
end

try
    X = load(Xname, 'X');
    X = X.X;
    assert(length(size(X)) == 2);
    assert(size(X, 2) == 3);
catch me,
    disp('Argument of this Matlab function must be a string corresponding to a .mat file (e.g. argument ''input.mat'').');
    disp('This .mat file should contain a variable X of size Nx3 describing the coordinates of the set of N points.');
    disp(' ');
    disp('Return values Vopt and Ropt will be stored in the file ''output.mat''.');
    disp('Vopt = Minimal OBB volume of the set of points (tolerance = 1e-9).');
    disp('Ropt = Rotation matrix associated to the optimal OBB.');
    disp(' ');
    disp('Bugs and errors can be reported on https://github.com/chagome/OptimalOBB/issues.');
    Vopt = [];
    Ropt = [];
    return;
end

try
    tol = 1e-9;

    % Computation of the geometry of the convex hull
    X = X - ones(size(X, 1), 1)*mean(X);
    normFactor = max(sqrt(sum(X.^2, 2)));
    X = X./normFactor;
    X = X(unique(convhulln(X)), :);
    hullFaces = convhulln(X);
    hullEdges = listEdges(hullFaces, X);
    hullNodes = nodeToFaces(hullFaces, X);
    hullNormals = zeros(size(hullFaces,1), 3);
    for k=1:size(hullNormals, 1)
        hullNormals(k,:) = -normalToFaces(hullFaces(k, :), X);
    end
    i0 = 2;
    j0 = 1;
    nEdges = size(hullEdges,1);

    % Computation of the minimal oriented bounding box
    Vopt = Inf;

    for i=i0:max(nEdges, i0),
        normal1 = hullNormals(hullEdges(i, [3 4]), :);
        real_mean_normal1 = [.5 .5]*normal1;
        real_e1 = hullEdges(i,5:7);
        for j=j0:min(i-1,j0+nEdges-1),
            normal2 = hullNormals(hullEdges(j,[3 4]),:);
            real_mean_normal2 = [.5 .5]*normal2;
            real_e2 = hullEdges(j,5:7);

            % Computation of the frame and the range of theta1
            [phi, R_x1, theta1_min, theta1_max, rc1] = computeTheta(real_e1, real_e2, tol, normal1, normal2, real_mean_normal1);
            switching = (abs(phi) <= tol*1e6) & (theta1_max - theta1_min) < 1e-4;
            switching(:) = 0;
            if abs(phi) <= tol*1e6,
                continue
            end

            for k=1:length(theta1_min),
                if switching(k),
                    R_x = R_x2;
                    mean_normal1 = real_mean_normal2;
                    mean_normal2 = real_mean_normal1;
                    e1 = real_e2;
                    e2 = real_e1;
                else,
                    R_x = R_x1;
                    mean_normal1 = real_mean_normal1;
                    mean_normal2 = real_mean_normal2;
                    e1 = real_e1;
                    e2 = real_e2;
                end
                x = R_x(1,:);
                z = R_x(3,:);

                % Declaration of some useful arrays
                p1 = zeros(1,4);
                p2 = zeros(1,4);
                theta1_p2 = (theta1_max(k)+2*tol)*ones(1,4);
                nodes_p2 = cell(1,4);
                nodeExt = cell(1,4);
                nextNodeExt = zeros(1,4);

                % Computation of the smallest bounding box oriented along the axis determined by the first theta1
                theta1 = theta1_min(k);
                R_n = computeR_n(theta1,phi,R_x,mean_normal2);%(1,:));
                X_n = X*R_n';
                [minX_n,nodeMin] = min(X_n);
                [maxX_n,nodeMax] = max(X_n);
                V = prod(maxX_n-minX_n);
                if V < Vopt,
                    % Update of the optimal volume and rotation
                    Ropt = R_n;
                    Vopt = V;
                end
                prevNodeExt = [nodeMin nodeMax(3)];
                for l=1:3,
                    nodeExt{l} = find(X_n(:,l) <= minX_n(l)+tol*1e2*max(abs(minX_n(l)),1));
                end
                nodeExt{4} = find(X_n(:,3) >= maxX_n(3)-tol*1e2*max(abs(maxX_n(3)),1));

                % Loop on all theta1 crossing an edge on the gaussian sphere
                while theta1 < theta1_max(k) + tol,

                    % Computation of the p2's for -n1, -n2, -n3 and +n3
                    for n=1:4,
                        if p2(n) == 0,
                            for m = 1:length(nodeExt{n}),
                                faces = hullNodes{nodeExt{n}(m)};
                                nFaces = length(faces);
                                if nFaces > 0,
                                    if xor(n==3,phi > 0),
                                        faces = fliplr(faces);
                                    end
                                    if n<=2,
                                        nIter = nFaces;
                                    else,
                                        nIter = nFaces+1;
                                    end
                                    if p1(n) > 0,
                                        l0 = find(faces == p1(n));
                                        if n <= 2
                                            l0 = l0 + 1;
                                        end
                                    else,
                                        l0 = 1;
                                    end
                                    if isempty(l0),
                                        l0 = 1;
                                    end
                                    for ll=1:nIter,
                                        l = max(1+mod(l0-1+ll-1,nFaces),1);
                                        face = hullFaces(faces(l),:);
                                        normals = hullNormals([faces(l);faces(1+mod(l,nFaces))],:);
                                        if xor(n==3,phi > 0)
                                            normals = hullNormals([faces(l);faces(1+mod(l-2,nFaces))],:);
                                        end
                                        idxNode = find(face == nodeExt{n}(m), 1);
                                        edge = [face(idxNode) face(1+mod(idxNode,3))];
                                        % [i j n ll nFaces abs((X(edge(2),:)-X(edge(1),:))*R_n(3,:)') edge]
                                        a = X(edge(2),:)-X(edge(1),:);
                                        if norm(a) < tol
                                            continue;
                                        end
                                        if n==1
                                            c = [a(2)*e1(3)-a(3)*e1(2) a(3)*e1(1)-a(1)*e1(3) a(1)*e1(2)-a(2)*e1(1)];
                                            if norm(c)/norm(a) > tol,
                                                sigma = -sign(c*R_n(1,:)');
                                                sigma(sigma == 0) = 1;
                                                intersection = sigma*c/norm(c);
                                            else,
                                                continue;
                                            end
                                            newN1 = -intersection;
                                        else
                                            if n==2
                                                c = [a(2)*e2(3)-a(3)*e2(2) a(3)*e2(1)-a(1)*e2(3) a(1)*e2(2)-a(2)*e2(1)];
                                                if norm(c)/norm(a) > tol,
                                                    sigma = -sign(c*R_n(2,:)');
                                                    sigma(sigma == 0) = 1;
                                                    intersection = sigma*c/norm(c);
                                                else,
                                                    continue;
                                                end
                                                c = intersection;
                                                c = [c(2)*e1(3)-c(3)*e1(2) c(3)*e1(1)-c(1)*e1(3) c(1)*e1(2)-c(2)*e1(1)];
                                                if norm(c)/norm(intersection) > tol,
                                                    sigma = sign(c*R_n(1,:)');
                                                    sigma(sigma == 0) = 1;
                                                    newN1 = sigma*c/norm(c);
                                                else,
                                                    continue;
                                                end
                                            else,
                                                edge_x = a*R_x';
                                                a = edge_x(1)*cos(phi)/2;
                                                b = -edge_x(2)*cos(phi)/2;
                                                c = edge_x(3)*sin(phi);

                                                psi = atan2(b,a);
                                                r = sqrt(a^2 + b^2);
                                                if abs(r) < max(abs(c-b),tol),
                                                    continue;
                                                end
                                                theta1inter = 0.5*([0 pi]+asin((c-b)/r)*[1 -1]-psi);
                                                theta1inter = mod([theta1inter theta1inter+pi],2*pi);
                                                if mean_normal1*x' >= 0,
                                                    theta1inter = theta1inter-pi;
                                                end
                                                idxTh1 = find(theta1inter > theta1);
                                                if isempty(idxTh1),
                                                    theta1inter = min(theta1inter);
                                                else,
                                                    theta1inter = min(theta1inter(idxTh1));
                                                end
                                                newR_n = computeR_n(theta1inter,phi,R_x,mean_normal2);
                                                newN1 = newR_n(1,:);
                                                intersection = (2*n-7)*newR_n(3,:);
                                            end
                                        end
                                        if n <= 2 || ll > 1 || abs((X(edge(2),:) - X(edge(1),:))*R_n(3,:)') > tol,
                                            a = [normals(1,:); intersection];
                                        else
                                            a = [(2*n-7)*R_n(3,:); intersection];
                                        end
                                        b = [intersection; normals(2,:)];
                                        c = [a(:,2).*b(:,3)-a(:,3).*b(:,2) a(:,3).*b(:,1)-a(:,1).*b(:,3) a(:,1).*b(:,2)-a(:,2).*b(:,1)];
                                        if norm(diff(a)) < tol,
                                            testBinary = 1;
                                        else
                                            testBinary = c(1,:)*c(2,:)' >= -tol;
                                        end
                                        if testBinary,
                                            a = R_n(1,:);
                                            b = newN1;
                                            c = [a(:,2).*b(:,3)-a(:,3).*b(:,2) a(:,3).*b(:,1)-a(:,1).*b(:,3) a(:,1).*b(:,2)-a(:,2).*b(:,1)];
                                            sinDtheta1 = c*z';
                                            cosDtheta1 = a*b';
                                            dtheta1 = real(atan2(sinDtheta1,cosDtheta1));
                                            if dtheta1 > tol,
                                                nextTheta1 = min(theta1+dtheta1,theta1_p2(n))-tol;
                                                midR_n = computeR_n(0.5*(theta1+nextTheta1),phi,R_x,mean_normal2);
                                                nextR_n = computeR_n(nextTheta1,phi,R_x,mean_normal2);
                                                delta2 = 1;
                                                if n<=3,
                                                    delta = -diff(X([prevNodeExt(n);edge(1)],:)*midR_n(n,:)');
                                                    if p2(n) > 0
                                                        delta2 = -diff(X([nextNodeExt(n);edge(2)],:)*nextR_n(n,:)');
                                                    end
                                                    delta3 = -diff(X([edge(1);edge(2)],:)*midR_n(n,:)');
                                                else,
                                                    delta = diff(X([prevNodeExt(n);edge(1)],:)*midR_n(3,:)');
                                                    if p2(n) > 0
                                                        delta2 = diff(X([nextNodeExt(n);edge(2)],:)*nextR_n(3,:)');
                                                    end
                                                    delta3 = diff(X([edge(1);edge(2)],:)*midR_n(3,:)');
                                                end
                                                if p2(n)==0 || p1(n) > 0 || (p1(n)==0 && p2(n)>0 && (theta1+dtheta1 < theta1_p2(n)-tol || (delta >=0 && delta2 > -tol && delta3 < tol)))
                                                    if theta1+dtheta1 < theta1_max(k) && (delta >=0 && delta2 > -tol && delta3 < tol),
                                                        prevNodeExt(n) = edge(1);
                                                        prevNodeExt_n  = prevNodeExt(n);
                                                        if xor(n==3,phi > 0)
                                                            p2(n) = faces(1+mod(l-2,nFaces));
                                                        else
                                                            p2(n) = faces(1+mod(l,nFaces));
                                                        end
                                                        theta1_p2(n) = theta1+dtheta1;
                                                        intersection = intersection*sign(intersection*mean(normals)');
                                                        nodes_p2{n} = findEquivalentExtrema(intersection,edge(2),hullFaces,hullNodes,X,tol*1e2,prevNodeExt_n);
                                                        nextNodeExt(n) = edge(2);
                                                        if nextNodeExt(n) == 0
                                                            assert(false);
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            if p2(n)==0,
                                X_n = X*R_n';
                                x_n = X*computeR_n(theta1_max(k),phi,R_x,mean_normal2)';
                                if n<=3,
                                    minX_n = min(X_n(:,n));
                                    [min_x,nodemin] = min(x_n(:,n));
                                    if any(find(X_n(:,n) < minX_n+tol*max(1e3*abs(minX_n),1)) == nodemin)
                                        prevNodeExt(n) = nodemin;
                                    end
                                else,
                                    maxX_n = max(X_n(:,3));
                                    [max_x,nodemax] = max(x_n(:,3));
                                    if any(find(X_n(:,3) > maxX_n-tol*max(1e3*abs(maxX_n),1)) == nodemax)
                                        prevNodeExt(4) = nodemax;
                                    end
                                end
                            end
                        end
                    end
                    p2(p2==0) = -1;
                    theta1_new = min([theta1_p2 theta1_max(k)]);
                    update = find(theta1_p2 <= theta1_new + tol);
                    p1(update) = p2(update);
                    p2(update) = 0;
                    theta1_p2(update) = theta1_max(k)+2*tol;

                    % Find the minimum volume within the range ]theta1,theta1_new]
                    X_Min = X(prevNodeExt(1:3),:);
                    if switching(k),
                        X_Max = [X(hullEdges([j i],1),:);X(prevNodeExt(4),:)];
                    else,
                        X_Max = [X(hullEdges([i j],1),:);X(prevNodeExt(4),:)];
                    end
                    c = X_Max-X_Min;
                    c1 = c(1,:) - (c(1,:)*e1')*e1;
                    gamma1 = -acos(-c1*x'/norm(c1));
                    c2 = c(2,:) - (c(2,:)*e2')*e2;
                    gamma2 = pi-real(asin([x(2)*c2(3)-x(3)*c2(2) x(3)*c2(1)-x(1)*c2(3) x(1)*c2(2)-x(2)*c2(1)]*e2'/norm(c2)));
                    c3 = c(3,:)*R_x';
                    num = sign(sin(phi))*[c3(3)*sin(phi) -c3(1)*cos(phi) c3(2)*cos(phi)+c3(3)*sin(phi)];
                    num = [num 0]*sin(gamma1)+[0 num]*cos(gamma1);
                    num = [num 0]*cos(gamma2)*sin(phi)+[0 num]*sin(gamma2);
                    deg = length(num)-1;
                    deriv = (deg:-1:1).*num(1:deg);
                    num = [deriv 0 0 0 0]*sin(phi)^2+[0 0 deriv 0 0]*(sin(phi)^2+1)+[0 0 0 0 deriv]-(4*sin(phi)^2*[num 0 0 0]+2*(sin(phi)^2+1)*[0 0 num 0]);
                    rootsNum = roots(num);
                    theta1Local = mod(atan(rootsNum(abs(imag(rootsNum))==0)),pi)-pi;
                    theta1Local = theta1Local(theta1Local >= theta1);
                    theta1 = [theta1Local(theta1Local < theta1_new);theta1_new];
                    for m=1:length(theta1)
                        R_n = computeR_n(theta1(m),phi,R_x,mean_normal2);
                        if m < length(theta1),
                            V = prod(diag(X_Max*R_n')-diag(X_Min*R_n'));
                        else,
                            X_Min = min(X([prevNodeExt(1:3) nextNodeExt(nextNodeExt(1:3) > 0)],:)*R_n');
                            X_Max = max(X([reshape(hullEdges([i j],1:2),1,4) prevNodeExt(4) nextNodeExt(4)*ones(nextNodeExt(4)>0,1)],:)*R_n');
                            V = prod(X_Max-X_Min);
                        end
                        if V < Vopt,
                            % Update of the optimal volume and rotation
                            Ropt = R_n;
                            Vopt = V;
                        end
                    end

                    for m=1:length(update),
                        n = update(m);
                        nodeExt{n} = nodes_p2{n};
                        prevNodeExt(n) = nextNodeExt(n);
                    end

                    theta1 = theta1_new+tol;
                end
            end
        end
    end
    Vopt = Vopt * normFactor^3;
    save('output.mat', 'Vopt', 'Ropt');
    fprintf('Vopt = %.10f\n', Vopt);
catch me,
    Vopt = [];
    Ropt = [];
    disp('An unexpected error has occured. Please report this on https://github.com/chagome/OptimalOBB/issues with the datafile.');
end



function c = cross(a,b)

c = a(:, [2 3 1]).*b(:, [3 1 2]) - a(:, [3 1 2]).*b(:, [2 3 1]);

function nodes = findEquivalentExtrema(normal,todo,allFaces,allNodes,X,tol,prev)

init = false(size(allNodes, 1), 1);
key = 1:size(allNodes, 1);
done = init;
done(todo) = 1;
keep = done;
done(prev) = 1;
limit = X(todo,:)*normal'-tol*max(abs(X(todo,:)*normal'),1);
while ~isempty(todo),
    tocheck = allFaces(cell2mat(allNodes(todo)'), :);
    processing = init;
    processing(tocheck) = 1;
    processing = key(processing & ~done);
    todo = processing(X(processing,:)*normal' >= limit);
    keep(todo) = 1;
    done(processing) = 1;
end
nodes = key(keep);

function [phi, R_x, theta1_min, theta1_max, rc] = computeTheta(e1, e2, tol, normal1, normal2, mean_normal1)

% Computation of the frame axes and angle
c = [e2(2)*e1(3)-e2(3)*e1(2) e2(3)*e1(1)-e2(1)*e1(3) e2(1)*e1(2)-e2(2)*e1(1)];
if norm(c) < tol*1e-1,
    c = [0 e1(3) -e1(2)];
    if norm(c) < tol*1e-1,
        c = [e1(2) -e1(1) 0];
    end
end
x = c/norm(c);
c = [e1(2)*x(3)-e1(3)*x(2) e1(3)*x(1)-e1(1)*x(3) e1(1)*x(2)-e1(2)*x(1)];
y = c/norm(c);
z = e1;
R_x = [x;y;z];
phi = real(asin([y(2)*e2(3)-y(3)*e2(2) y(3)*e2(1)-y(1)*e2(3) y(1)*e2(2)-y(2)*e2(1)]*x'));

% Computation of the bound on theta1
rc = 0;
if abs(phi) > tol*1e6,
    n1 = [normal1;cross(normal2,[e1;e1])];
    sigma = n1(3:4,:)*mean_normal1';
    sigma(abs(sigma) < tol) = 1;
    sigma = sign(sigma)./[norm(n1(3,:));norm(n1(4,:))];
    n1(3:4,:) = n1(3:4,:).*(sigma*[1 1 1]);
    theta1 = atan2(cross([x;x;x;x],n1)*e1',n1*x');
else,
    c = cross(normal2,[z;z]);
    if mean_normal1*y' > 0,
        sigma = 1;
    else,
        if mean_normal1*x' < 0,
            sigma = 3;
        else,
            sigma = -1;
        end
    end
    if c(1,:)*c(2,:)' < -tol*1e2,
        theta1 = atan2(cross([x;x],normal1)*e1',normal1*x');
        c = cross(normal1,[y;y]);
        if c(1,:)*c(2,:)' < -tol,
            rc = sigma;
        else,
            theta1 = [atan2(cross([x;x],normal1)*e1',normal1*x');0;2*pi];
        end
    else,
        theta1 = [atan2(cross([x;x],normal1)*e1',normal1*x');sigma*pi/2;sigma*pi/2];
    end
end
if mean_normal1*x' < 0,
    theta1(theta1<0) = theta1(theta1<0) + 2*pi;
end
if rc,
    theta1_min = [theta1(1) rc*pi/2];
    theta1_max = [rc*pi/2 theta1(2)];
else,
    if phi < 0
    theta1(3:4) = theta1([4 3]);
    end
    if theta1(4) >= theta1(3) - tol
    if theta1(2) - theta1(1) < -tol
        assert(false);
    end
    theta1_min = max(theta1([1 3]));
    theta1_max = min(theta1([2 4]));

    else,
    theta1_min = max(theta1([1 1;3 3]) - [0 0;pi 0]);
    theta1_max = min(theta1([2 2;4 4]) + [0 0;0 pi]);
    end
end

gud = theta1_min < theta1_max-tol*1e3;
theta1_min = theta1_min(gud)+tol;
theta1_max = theta1_max(gud)-tol*2e1;
