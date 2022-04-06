function [X,res_d,res_p]= SRSSWSU_1(Y,A,lambda_l1,mu)


norm_y = sqrt(mean(Y(:).^2));
Y = Y./norm_y;
A = A./norm_y;
MaxIter_1 = 50;
MaxIter_2 = 5;
epsilon = 1e-5;
[L, N] = size(Y);
m = size(A, 2);
Finv = inv(A'*A+2*eye(m));
U = Finv*A'*Y;

%Initialization of auxiliary matrices V1,V2,V3,V4
V1 = A*U; 
V2 = U;
V3 = U;
%Initialization of Lagranfe Multipliers
D1 = V1*0;
D2 = V2*0;
D3 = V3*0;

%Initialization of iteration number
i = 1;
K = 1;
%primal residual 
res_p = inf;
%dual residual
res_d = inf;

%error tolerance
toll = sqrt((3*m + L)/2*N/2)*epsilon;

while (i <= MaxIter_1)
Temp = sqrt(sum((V2-D2).^2,2));
b=1./Temp;   

w1=repmat(b,1,size(V2,2));

%while (K <= MIter_2) && ((abs(res_p) > toll) || (abs(res_d) > toll))
while (K <= MaxIter_2)
    if mod(K, 5) == 1
        V10 = V1;
        V20 = V2;
        V30 = V3;
    end
    
    Wprev = U;
    
    U = Finv*(A'*(V1 + D1) + (V2 + D2) + (V3 + D3));
    V1 = 1/(1+mu)*(Y + mu*(A*U - D1));
    
    
    V2 = soft(U-D2,lambda_l1/mu.*w1);
    %positivity             
    V3 = max(U - D3, 0);
    
    %update D
    D1 = D1 - A*U + V1;
    D2 = D2 - U + V2;
    D3 = D3 - U + V3;
     
    if mod(K, 5) == 1
        %primal residual
        res_p = norm([V1; V2; V3] - [A*U; U; U], 'fro');
        %dual residual
        res_d = norm([V1; V2; V3] - [V10; V20; V30], 'fro');
       
    end
   % error(i) = norm(W - Wtrue,'fro');
    K = K + 1;
end

X = U;
i=i+1;
K=1;
end


