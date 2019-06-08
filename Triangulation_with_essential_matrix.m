close all
clear all
addpath('functions/')
addpath('images/')
addpath('data/')

% Loading data
load('ex5data.mat')

% Normalized the points
Nx1 = pflat(K\x{1});
Nx2 = pflat(K\x{2});

% Computing the Essential matrix using the using 3 point-matches
% and RANSAC algorithm to select the best solution.
[E,inliers,U,V] = RANSAC(Nx1,Nx2,K,600);

% Computing camera matrices
W = [0 -1 0; 1 0 0; 0 0 1];
P1 = [eye(3) zeros(3,1)];
U3 = U(:,3);

P2{1} = [U*W*V' U3];
P2{2} = [U*W*V' -U3];
P2{3} = [U*W'*V' U3];
P2{4} = [U*W'*V' -U3];

N = [];
for i = [1 2 3 4]
    % Triangulating by using all four camera matrices P2
    X{i} = triangulation(Nx1(:,inliers),Nx2(:,inliers),P1,P2{i});
    
    % Projecting the triangulated scene point X on camera P1
    xp{i} = P1*X{i};
    
    % Projecting the triangulated scene point X on camera P2
    x2p{i} = P2{i}*X{i};
    
    % Counting all the points in front of camera pairs P1 and P2
    N = [N sum(x2p{i}(3,:)>0)+sum(xp{i}(3,:)>0)];
end
% Extracting the index for the correct camera
% which should have most scene points in front
% of it.
[value,index]=max(N);

% Solution: Scene point X and Camera matrix P2
X = X{index};
P2 = P2{index};

% Plotting the triangulated scene points
plot3(X(1,:),X(2,:),X(3,:),'.');

%% Functions

% Computes three matrices E1, E2 and E3
function E123 = computeEs(x1,x2)
    X1=x1(1,:)';
    Y1=x1(2,:)';
    Z1=x1(3,:)';
    X2=x2(1,:)';
    Y2=x2(2,:)';
    Z2=x2(3,:)';
    
    % Creating the M matrix used to compute E1, E2 and E3
    M = [X2.*X1-Y2.*Y1 X2.*Y1-Y2.*X1 X2.*Z1  Y2.*Z1 Z2.*X1 Z2.*Y1];
    
    % Computing the non-trivial solution to nullspace of M
    e = null(M);
    
    % Reshaping the elements of the solutions into matrices E1, E2 and E3.
    E123 = {};
    for i = [1 2 3]
        E123{i} = [e(1,i)  e(2,i) e(3,i)
             e(2,i) -e(1,i) e(4,i)
             e(5,i)  e(6,i)     0];
    end
end

function [bestE,bestInliers,bestU,bestV] = RANSAC(y,x,K,iterations)

    bestAmountOfInliers = 0;
    for i = 1:iterations
        
        % Choosing 3 random matching pair of points
        perm = randperm(size(x ,2));
        xrand = x(:,perm(1:3));
        yrand = y(:,perm(1:3));
        
        % Computing essential matrices E1, E2 and E3 using 3 number of correspondences
        E123 = computeEs(yrand,xrand);
        E = solve_polynomial_sys(E123);
        
        % Count inliers and return the corresponding one with the
        % highest amount
        for i = 1:length(E)
            
            [U, S, V] = svd(E{i});
            %Check minimum singular value and Mv are both small

            % Constructing Essential matrix
            if det(U*V') >0
                E{i} = U*diag([1 1 0])*V';
            else
                V = -V;
                E{i} = U*diag([1 1 0])*V';
            end
            
            % Count inliers
            [amountOfInliers, inliers] = countInliers(y,x,E{i},K);
        
            % Saving the camera matrix P with most inliers
            if amountOfInliers > bestAmountOfInliers
                bestE = E{i};
                bestInliers = inliers;
                bestAmountOfInliers = amountOfInliers;
                bestU = U;
                bestV = V;
                disp(strcat('Amount of inliers-  ',num2str(bestAmountOfInliers)))
            end
        end
    end
end

function [amount, inliers] = countInliers(y,x,Essential,K)
        
    % Fundamental matrix
    F = (K'\Essential)*(K\eye(3));
    F = F./F(3,3);

    % Epipolar line from x
    l = F*K*x;

    % Calculate distace from point and line, distances larger than 10 
    % count as outliers
    amount = sum((abs((sum(l.*(K*y),1))./sqrt(l(1,:).^2+l(2,:).^2)))<10);
    inliers = (abs((sum(l.*(K*y),1))./sqrt(l(1,:).^2+l(2,:).^2)))<10;
end

function X = triangulation(x1,x2,P1,P2) 
    xm1 = mean(x1(1:2,:),2);
    xm2 = mean(x2(1:2,:),2);
    s1 = std(x1(1:2,:),0,2);
    s2 = std(x2(1:2,:),0,2);

    % Normalization matrices
    N1 = [1/mean(s1) 0          -(1/mean(s1))*xm1(1)
          0          1/mean(s1) -(1/mean(s1))*xm1(2)
          0          0          1];
    N2 = [1/mean(s2) 0          -(1/mean(s2))*xm2(1)
          0          1/mean(s2) -(1/mean(s2))*xm2(2)
          0          0          1];

    X = [];
    for i = 1:length(x1(1,:))
        M = [N1*P1;N2*P2];
        M = [M,[-N1*x1(:,i) zeros(3,1);zeros(3,1) -N2*x2(:,i)]];
        [~,~,V] = svd(M);
        X = [X, V(1:4,end)];
    end
    
    X=pflat(X);
end
