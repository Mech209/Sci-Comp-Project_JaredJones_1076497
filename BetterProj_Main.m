clear;clc;close all
%Domain
ax = -pi();
ay = -pi();
bx = pi();
by = pi();
tol = 1e-6; %degree of tolerance
it = 50; %Max Number of Iterations to complete
%Delta
N = 3;
h = (bx-ax)/(N+1);
lambda = 1;
x = ax:h:bx;
x = x';
y = ay:h:by;
y = y';
D = (-4 + lambda*h^2);

%Building U matrix with Boundary Conditions and Given F 
u = zeros(N+2);

for k = 1:N+2
    u(k,N+2) = (bx - x(k))^2 * cos(pi()*x(k)/(bx)); %BC u(x,y=by) = fb(x)
    u(k,1) = x(k)*(bx - x(k))^2; %BC u(x,y=ay) = gb(x)
    u(1,k) = (bx - ax)^2 * cos(pi()*ax/(bx)) + ((y(k) - ay)/(by - ay))*(ax*(bx -ax)^2 - (bx - ax)^2 * cos(pi()*ax/(bx)));
end
F = zeros();
for j = 1:N+2
   for k = 1:N+2
      F(k,j) = (sin(pi().*((x(k)-ax)/(bx-ax))).*cos((pi()/2).*(2.*((y(j)-ay)/(by-ay)) +1))); 
   end
end

%RHS Vector, F
F = F(2:N+2,2:N+1);
F(:,1) = F(:,1) - u(2:N+2,1);
F(1,:) = F(1,:) - u(1,2:N+1);
F(:,N) = F(:,N) - u(2:N+2,N+2);
F = reshape(F,(N+1)*N,1);

%Uknown U Matrix
U = u(2:N+2,2:N+1);
U = reshape(U,(N+1)*N,1);

%creating coefficient matrix
K_sup = eye(N+1);
K_Main = D*eye(N+1) + diag(ones(N,1),1) + diag(ones(N,1),-1); 
K_Main(N+1,N) = 2;

K = sparse((N+1)*N,(N+1)*N);
for ii = 1:N
    K(1+(ii-1)*(N+1):(N+1)+(ii-1)*(N+1),1+(ii-1)*(N+1):(N+1)+(ii-1)*(N+1))= K_Main;
end

for jj = 1:N-1
    K(jj*(N+1)+1:jj*(N+1)+(N+1),1+(jj-1)*(N+1):(N+1)+(jj-1)*(N+1))=K_sup;
end

for kk = 1:N-1
    K(1+(kk-1)*(N+1):(N+1)+(kk-1)*(N+1),kk*(N+1)+1:kk*(N+1)+(N+1))=K_sup;
end

A = GaussSeidel(K,F,tol,it)