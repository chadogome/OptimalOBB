function [Ropt,Vopt] = caliper3D(X,hullFaces,hullEdges,hullNodes,hullNormals,i,j)

% @pre    : X(i,j) is the j^th coordinate of the node with index i
% @pre    : hull***** is defined in oRourke.m
% @return : R_opt(i,:) is the i^th axis for the frame aligned with the minimum volume oriented bounding box
% @return : V_opt is the minimal volume over all the oriented bounding box
%
% MELCHIOR Samuel - UCL, 2008-09.

tol = 1e-10

% Computation of the minimal oriented bounding box
Vopt = Inf;

normal1 = hullNormals(hullEdges(i,[3 4]),:);
e1 = hullEdges(i,[5:7]);
normal2 = hullNormals(hullEdges(j,[3 4]),:);
e2 = hullEdges(j,[5:7]);

% Computation of the frame axes and angle
c = [e2(2)*e1(3)-e2(3)*e1(2) e2(3)*e1(1)-e2(1)*e1(3) e2(1)*e1(2)-e2(2)*e1(1)];
x = c/norm(c);
c = [e1(2)*x(3)-e1(3)*x(2) e1(3)*x(1)-e1(1)*x(3) e1(1)*x(2)-e1(2)*x(1)];
y = c/norm(c);
z = e1;
R_x = [x;y;z];
phi = asin(cross(y,e2)*x');

% Computation of the bound on theta1
sigma = sign(normal1*y');
sigma(sigma==0)=1;
n1 = -diag(sigma)*normal1;
theta1 = -real(acos(n1*x'));      % = asin(cross(x,n1)*e1');
sigma = sign(normal2*x');
sigma(sigma==0)=1;
n2 = diag(sigma)*normal2;
theta2 = real(asin(cross([1;1]*x,n2)*e2'));
theta1(:,2) = theta2toTheta1(theta2,phi);
if phi < 0
    theta1(:,2) = flipud(theta1(:,2));
end
if (theta1(1,1) - theta1(2,1))*(theta1(1,2) - theta1(2,2)) >= 0
    theta1_min = max(theta1(1,:));
    theta1_max = min(theta1(2,:));
else
    if theta1(1,1) > theta1(2,1)
        theta1_min = max([theta1(1,:);theta1(1,:)] - [pi 0;0 0]);
        theta1_max = min([theta1(2,:);theta1(2,:)] + [pi 0;0 0]);
    else
        theta1_min = max([theta1(1,:);theta1(1,:)] - [0 pi;0 0]);
        theta1_max = min([theta1(2,:);theta1(2,:)] + [0 pi;0 0]);
    end
end
gud = theta1_min <= theta1_max;
theta1_min = mod(theta1_min(gud),pi)-pi;
theta1_max = mod(theta1_max(gud),pi)-pi;

% Loop on all theta1 crossing an edge on the gaussian sphere
for k=1:length(theta1_min)
    theta1 = theta1_min(k);
    p1 = zeros(1,4);
    p2 = zeros(1,4);
    nNew = [];
    nodeExt = [];

    while theta1 <= theta1_max(k)
        % Computation of the smallest bounding box oriented along the axis determined by theta1
        R_n = computeR_n(theta1,phi,R_x);
        X_n = X*R_n';
%         X_n(nodeExt(nNew),min(nNew,3)) = mean(X_n(:,min(nNew,3)));    COULD I COMMENT IT ????
        [X_nMin,nodeMin] = min(X_n);
        [X_nMax,nodeMax] = max(X_n);
        V = prod(X_nMax-X_nMin);
        if V < Vopt
            % Update of the optimal volume and rotation
            Ropt = R_n;
            Vopt = V;
        end

        % Computation of the p2's for -n1, -n2, -n3 and +n3
        if R_n(1,:)*normal1(1,:)' > 0;
            nodeExt = [nodeMin nodeMax(3)];
        else
            nodeExt = [nodeMax nodeMin(3)];
        end
        for n=1:4
            if p2(n) == 0
                faces = hullNodes{nodeExt(n)};
                nFaces = length(faces);
                if n > 2
                    l0 = p1(n);
                else
                    l0 = p1(n) + 1;
                end
                for ll=1:nFaces
                    l = max(1+mod(l0-1+ll-1,nFaces),1);
                    face = hullFaces(faces(l),:);
                    normals = hullNormals([faces(l);faces(1+mod(l-2,nFaces))],:);

                    idxNode = find(face == nodeExt(n));
                    edge = [face(idxNode) face(1+mod(idxNode,3))];
                    if n==1
                        a = X(edge(2),:)-X(edge(1),:);
                        c = [a(2)*e1(3)-a(3)*e1(2) a(3)*e1(1)-a(1)*e1(3) a(1)*e1(2)-a(2)*e1(1)];
                        intersection = c/norm(c);
                        newN1 = -intersection;
                    else
                        if n==2
                            a = X(edge(2),:)-X(edge(1),:);
                            c = [a(2)*e2(3)-a(3)*e2(2) a(3)*e2(1)-a(1)*e2(3) a(1)*e2(2)-a(2)*e2(1)];
                            intersection = c/norm(c);
                            c = [c(2)*e1(3)-c(3)*e1(2) c(3)*e1(1)-c(1)*e1(3) c(1)*e1(2)-c(2)*e1(1)];
                            newN1 = c/norm(c);
                            sigma = sign(newN1*R_n(1,:)');
                            sigma(sigma == 0) = 1;
                            newN1 = newN1*sigma;
                        else
                            edge_x = (X(edge(2),:)-X(edge(1),:))*R_x';
                            a = edge_x(1)*cos(phi)/2;
                            b = -edge_x(2)*cos(phi)/2;
                            c = edge_x(3)*sin(phi);
                            psi = atan2(b,a);
                            r = sqrt(a^2 + b^2);
                            if abs(r) < abs(c+b)
                                continue;       % there is no intersection
                            end
                            theta1inter = mod(0.5*([0 pi]+asin(-(c+b)/r)*[1 -1]-psi),pi)-pi;
                            idxTh1 = find(theta1inter > theta1);
                            if length(idxTh1) == 0
                                theta1inter = min(theta1inter);
                            else
                                theta1inter = min(theta1inter(idxTh1));
                            end
                            newR_n = computeR_n(theta1inter,phi,R_x);
                            newN1 = newR_n(1,:);
                            intersection = (2*n-7)*newR_n(3,:);
                        end
                    end
                    a = [normals(1,:); intersection];
                    b = [intersection; normals(2,:)];
                    c = [a(:,2).*b(:,3)-a(:,3).*b(:,2) a(:,3).*b(:,1)-a(:,1).*b(:,3) a(:,1).*b(:,2)-a(:,2).*b(:,1)];
                    testBinary = c(1,:)*c(2,:)' >= 0;
                    if testBinary
                        a = R_n(1,:);
                        b = newN1;
                        c = [a(:,2).*b(:,3)-a(:,3).*b(:,2) a(:,3).*b(:,1)-a(:,1).*b(:,3) a(:,1).*b(:,2)-a(:,2).*b(:,1)];
                        sinDtheta1 = c*z';
                        dtheta1 = asin(sinDtheta1); %cross(R_n(1,:),newN1)*z');

                        if dtheta1 > tol
                            if p2(n)==0 || p1(n) > 0 || (p1(n)==0 && p2(n)>0 && theta1+dtheta1 < theta1_p2(n))
                                theta1_p2(n) = theta1+dtheta1;
                                p2(n) = l;
                            end
                            if p1(n) > 0
                                break;
                            end
                        else
                            p1(n) = l;
                        end
                    end
                end
            end
        end
        p2(p2==0) = -1;
        [theta1_new,nNew] = min(theta1_p2);
        update = find(theta1_p2 <= theta1_new + tol);
        p1(update) = p2(update);
        p2(update) = 0;
        theta1_p2(update) = Inf;

        % Find the minimum volume within the range [theta1,theta1_new]
        c = X(nodeMax,:)-X(nodeMin(1:3),:);
        c1 = c(1,:) - (c(1,:)*e1')*e1;
        gamma1 = -acos(-c1*x'/norm(c1));      % = asin(cross(x,c1)*e1'/norm(c1));
        c2 = c(2,:) - (c(2,:)*e2')*e2;
        gamma2 = pi-asin(cross(x,c2)*e2'/norm(c2));
        c3 = c(3,:)*R_x';
        num = sign(sin(phi))*[c3(3)*sin(phi) -c3(1)*cos(phi) c3(2)*cos(phi)+c3(3)*sin(phi)];    % d3*R/cos^2 [R^2/cos^2 = (1+x^2*sin(phi)^2); cos^2 = 1/(1 + x^2)]
        num = [num 0]*sin(gamma1)+[0 num]*cos(gamma1);                                          % d3*d1*R/(cos^3*norm(c1))
        num = [num 0]*cos(gamma2)*sin(phi)+[0 num]*sin(gamma2);                                 % d3*d1*d2*R^2/(cos^4*norm(c1)*norm(c2))
        deg = length(num)-1;
        deriv = (deg:-1:1).*num(1:deg);
        num = [deriv 0 0 0 0]*sin(phi)^2+[0 0 deriv 0 0]*(sin(phi)^2+1)+[0 0 0 0 deriv]-(4*sin(phi)^2*[num 0 0 0]+2*(sin(phi)^2+1)*[0 0 num 0]);
        rootsNum = roots(num);
        theta1Local = mod(atan(rootsNum(abs(imag(rootsNum))==0)),pi)-pi;
        theta1Local = theta1Local(theta1Local >= theta1_min(k));
        theta1Local = theta1Local(theta1Local <= theta1_max(k));
        for n=1:length(theta1Local)
            R_n = computeR_n(theta1Local(n),phi,R_x);
            X_n = X*R_n';
            X_nMin = min(X_n);
            X_nMax = max(X_n);
            V = prod(X_nMax-X_nMin);
            if V < Vopt
                % Update of the optimal volume and rotation
                Ropt = R_n
                Vopt = V
            end
        end

        % Prepare next iteration
        theta1_new-theta1;
        theta1 = theta1_new;
    end
end