function Esol = solve_polynomial_sys(E)
%Finds Essential matrices such Esol{i} = x*E{1}+y*E{2}+E{3}, where E is a
%cell containing 3 3x3 matrices.

monMat = NaN(4,4);
monMat(4,1) = 1; %x^3
monMat(3,2) = 2; %x^2*y
monMat(2,3) = 3; %x*y^2
monMat(1,4) = 4; %y^3,
monMat(3,1) = 5; %x^2,
monMat(2,2) = 6; %x*y,
monMat(1,3) = 7; %y^2,
monMat(2,1) = 8; %x,
monMat(1,2) = 9; %y
monMat(1,1) = 10; %1
monMap = @(a,b)monMat(a+1,b+1);


% The constraint: 2*E*E'*E-trace(E*E')*E = 0
coeffs = zeros(6,10);
for i = 1:3
    for j = 1:3
        for k = 1:3
            col = monMap(sum([i,j,k]==1),sum([i,j,k]==2));
            
            new_coeffs = 2*E{i}*E{j}'*E{k}-trace(E{i}*E{j}')*E{k};
            new_coeffs = new_coeffs(2:3,:); %2nd and third row rank 6 basis
            coeffs(:,col) = coeffs(:,col) + new_coeffs(:);
        end
    end
end

% The determinant constraint
% e11*e22*e33+e12*e23*e31+e13*e21*e32-e11*e23*e32-e12*e21*e33-e13*e22*e31=0
det_coeffs = zeros(1,10);
for i = 1:3
    for j = 1:3
        for k = 1:3
            col = monMap(sum([i,j,k]==1),sum([i,j,k]==2));
            
            det_coeffs(:,col) = det_coeffs(:,col)+...
            E{i}(1,1)*E{j}(2,2)*E{k}(3,3) ... %e11*e22*e33
            + E{i}(1,2)*E{j}(2,3)*E{k}(3,1) ... %e12*e23*e31
            + E{i}(1,3)*E{j}(2,1)*E{k}(3,2) ... %e13*e21*e32
            - E{i}(1,1)*E{j}(2,3)*E{k}(3,2) ... %e11*e23*e32
            - E{i}(1,2)*E{j}(2,1)*E{k}(3,3) ... %e12*e21*e33
            - E{i}(1,3)*E{j}(2,2)*E{k}(3,1);    %e13*e22*e31
        end
    end
end

%coeffs = [coeffs; det_coeffs]; %Determinant constraint unused?

reduced_coeffs = rref(coeffs);

%construct action matrix multiplication with x
Mx = [-reduced_coeffs([3 5 6],7:end); 0 1 0 0]; 

%Extract real solutions
[V,D] = eig(Mx);
V = V./repmat(V(end,:),[4 1]);

Esol = cell(0);
for i = 1:4
    if isreal(D(i,i))
        x = V(2,i);
        y = V(3,i);
        
        Esol{end+1} = E{1}*x+E{2}*y+E{3};
    end
end
