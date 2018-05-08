clear;clc;close all
%Domain
ax = -pi();
ay = -pi();
bx = pi();
by = pi();

%Delta
N = 30;
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

%creating coefficient matrix
K_sup = eye(N+1);
K_Main = D*eye(N+1) + diag(ones(N,1),1) + diag(ones(N,1),-1); 
K_Main(N+1,N) = 2;

K = zeros((N+1)*N,(N+1)*N);
for ii = 1:N
    K(1+(ii-1)*(N+1):(N+1)+(ii-1)*(N+1),1+(ii-1)*(N+1):(N+1)+(ii-1)*(N+1))= K_Main;
end

for jj = 1:N-1
    K(jj*(N+1)+1:jj*(N+1)+(N+1),1+(jj-1)*(N+1):(N+1)+(jj-1)*(N+1))=K_sup;
end

for kk = 1:N-1
    K(1+(kk-1)*(N+1):(N+1)+(kk-1)*(N+1),kk*(N+1)+1:kk*(N+1)+(N+1))=K_sup;
end

%Begin Gaussian Elimination sequence
%K is coefficient Matrix
%F is RHS vector
[m,n] = size(K);
nb = n + 1;
   Aug = [K, F];
   %Forward Elimination
   for k = 1:n-1
       for i = k+1:n
           factor = Aug(i,k)/Aug(k,k);
           Aug(i,k:nb) = Aug(i,k:nb) - factor*Aug(k,k:nb);
       end
   end
   %Back Substitution 
   w = zeros(n,1);
   w(n) = Aug(n,nb)/Aug(n,n);
   for i = n-1:-1:1
       w(i) = (Aug(i,nb) - Aug(i,i+1:n)*w(i+1:n))/Aug(i,i);
   end
w = reshape(w,N+1,N);
u(2:N+2,2:N+1) = w;

surf(y,x,u)