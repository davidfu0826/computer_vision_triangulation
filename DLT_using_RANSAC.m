clear all
close all
addpath('functions/')
addpath('images/')
addpath('data/')

% Loading data
im1 = imread('prison1.jpg');
im2 = imread('prison2.jpg');
load('ex4data.mat');

% Points in front of camera???

% Determines which of the points are visible in image 1 and 2
visible1 = isfinite(x{1}(1,:));
visible2 = isfinite(x{2}(1,:));

% DLT using RANSAC
amountOfIterations = 40;
P1 = RANSAC(X(:,visible1),x{1}(:,visible1),false,amountOfIterations);
P2 = RANSAC(X(:,visible2),x{2}(:,visible2),false,amountOfIterations);

% RK-factorization
[K1,R1,t1] = RQ(P1);
[K2,R2,t2] = RQ(P2);

% Camera centers and principal axes
C1 = -R1'*t1;
C2 = -R2'*t2;
v1 = 10*R1(:,3)/norm(R1(:,3));
v2 = 10*R2(:,3)/norm(R2(:,3));

%% Plot
hold on
xlim([-50,30])
ylim([-10,15])
zlim([0,40])
plot3(X(1,:),X(2,:),X(3,:), '.', 'MarkerSize', 2)
plot3(C1(1),C1(2),C1(3),'g.', 'MarkerSize', 30)
plot3(C2(1),C2(2),C2(3),'r.', 'MarkerSize', 30)
quiver3(C1(1),C1(2),C1(3),v1(1),v1(2),v1(3),'g', 'MarkerSize', 20, 'LineWidth', 1.5, 'MaxHeadSize', 8)
quiver3(C2(1),C2(2),C2(3),v2(1),v2(2),v2(3),'r', 'MarkerSize', 20, 'LineWidth', 1.5, 'MaxHeadSize', 8)

%% Functions

% Creating M-matrix used for Direct Linear Transformation (DLT)
function M = createM(X,x)
    l=length(X(1,:));
    M = [];
    for i = 1:l
        M = [M;[X(:,i)' zeros(1,4) zeros(1,4)
        zeros(1,4) X(:,i)' zeros(1,4)
        zeros(1,4) zeros(1,4) X(:,i)']];
    end
    for i = 1:l
        M = [M,[zeros(3*(i-1),1);-x(:,i);zeros(3*(l-i),1)]];
    end
end

% Performs Direct Linear Transformation (DLT)
function P = DLT(X,x)
    
    % Setting up the linear homogeneous system Mv=0
    M=createM(X,x); 

    %Computing the singular value decomposition
    [U,S,V] = svd(M);

    %Extracting solution v, from last column of V
    v=V(:,end);

    %Extracting the camera matrices from the solutions
    P = [v(1:4)' 
        v(5:8)'
        v(9:12)'];
end

% Using RANSAC algorithm, picking random set of points each iterations
% and computes the amount of inliers. The camera matrix P with largest
% amount of inliers is selected.
function P = RANSAC(X,x,normalize,iterations)

    % Normalization
    if(normalize)
        xm = mean(x(1:2,:),2);
        s = std(x(1:2,:),0,2);

        % Normalization matrices
        N = [1/mean(s)  0          -(1/mean(s))*xm(1)
              0          1/mean(s)  -(1/mean(s))*xm(2)
              0          0          1];
    else
        N = eye(3);
    end
    
    Nx = N*x;

    bestAmountOfInliers = 0;
    for i = 1:iterations
        
        % Choosing 6 random matching pair of points
        perm = randperm(size(Nx ,2));
        Nxrand = Nx(:,perm(1:6));
        Xrand = X(:,perm(1:6));
        
        % Computing camera matrix P
        P = N\DLT(Xrand,Nxrand);
        
        % Count inliers
        amountOfInliers = countInliers(X,x,P);
        
        % Saving the camera matrix P with most inliers
        if amountOfInliers > bestAmountOfInliers
            bestCamera = P;
            bestAmountOfInliers = amountOfInliers;
        end
    end
    bestAmountOfInliers
    
    % Checking that the majority of scene points X are in front of camera P
    P = checkSign(bestCamera,X);
end

% Amount of inliers estimation by calculating the distance
% between the projected points with the original
function amount = countInliers(X,x,P)
    amount = 0;
    for i = 1:length(x(1,:))
        
        % If the distance is less than 5 pixels then
        % the pair is an inlier
        if norm(pflat(P*X(:,i)) - x(:,i)) < 5
            amount = amount + 1;
        end
    end
end

function newP=checkSign(P,X)
    px=P*X;
    % A negative third coordinate equates to that point being behind
    % camera P. Thus we sum the total amount of points behind and in
    % front.
    totNeg = sum(px(3,:) < 0);
    totPos = sum(px(3,:) > 0);
    
    % If most scene points are behind the camera, change the sign
    % so that most of the scene points are in front of the camera
    if totNeg > totPos
        newP = -P;
    else
        newP = P;
    end
end

% RQ-factorization based on the fact that the K-matrix is upper triangular
% and that the rotational matrix R is orthogonal. Taking advantage of
% orthogonality, we can use scalar product and norm to compute matrices K
% and R. Thus computing the inner parameters of P
function [K,R,t] = RQ(P)
    f = norm(P(3,1:3));
    R3 = P(3,1:3)/norm(P(3,1:3));
    e = dot(R3,P(2,1:3));
    dR2 = P(2,1:3)-e*R3;
    d = norm(dR2);
    R2 = dR2/d;
    c = dot(R3,P(1,1:3));
    b = dot(R2,P(1,1:3));
    aR1 = P(1,1:3)-b*R2-c*R3;
    a = norm(aR1);
    R1 = aR1/a;
    
    K = [a b c
        0 d e
        0 0 f]./f;

    normalizedP = K\P;
    normalizedP = normalizedP./norm(normalizedP(:,1));
    
    R = normalizedP(:,1:3);
    t = normalizedP(:,4);
end
