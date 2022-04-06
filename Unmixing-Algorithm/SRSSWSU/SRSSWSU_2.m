function [X]= SRSSWSU_2(Y,A,lambda_l1,beta,mu,Xd)


norm_y = sqrt(mean(Y(:).^2));
Y = Y./norm_y;
A = A./norm_y;
MaxIter = 200;
epsilon = 1e-5;
[L, N] = size(Y);
m = size(A, 2);
XD = Xd;
Finv = inv(A'*A +  3*eye(m));
U = Finv*A'*Y;

%Initialization of auxiliary matrices V1,V2,V3,V4
V1 = A*U; 
V2 = U;
V3 = U;
V4 = U;
%Initialization of Lagranfe Multipliers
D1 = V1*0;
D2 = V2*0;
D3 = V3*0;
D4 = V4*0;

%Initialization of iteration number
K = 1;
%primal residual 
res_p = inf;
%dual residual
res_d = inf;

%error tolerance
toll = sqrt((3*m + L)/2*N/2)*epsilon;



%---------------------------------------------
%  ADMM iterations
%---------------------------------------------
%while (K <= MIter) && ((abs(res_p) > toll) || (abs(res_d) > toll))
while (K <= MaxIter)
    disp(K)
    if mod(K, 10) == 1
        V10 = V1;
        V20 = V2;
        V30 = V3;
        V40 = V4;
    end
    
    Wprev = U;
    
    U = Finv*(A'*(V1 + D1) + (V2 + D2) + (V3 + D3) + (V4 + D4));
    V1 = 1/(1+mu)*(Y + mu*(A*U - D1));
    
    
    V2 = soft(U-D2,lambda_l1/mu);
    V3 = (1./(beta+mu).*(beta.*XD+mu*(U-D3)));
    %positivity             
    V4 = max(U - D4, 0);
    
    %update D
    D1 = D1 - A*U + V1;
    D2 = D2 - U + V2;
    D3 = D3 - U + V3;
    D4 = D4 - U + V4;
    
    if mod(K, 10) == 1
        %primal residual
        res_p = norm([V1; V2; V3; V4] - [A*U; U; U; U], 'fro');
        %dual residual
        res_d = norm([V1; V2; V3; V4] - [V10; V20; V30; V40], 'fro');
   % error(i) = norm(W - Wtrue,'fro');
   
   end
K = K + 1;
X = U;

end
