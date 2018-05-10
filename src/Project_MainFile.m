clear;clc;close all
%Domain
ax = -pi();
ay = -pi();
bx = pi();
by = pi();

%Delta
N = 75;
h = (bx-ax)/(N+1);
x = ax:h:bx;
y = ay:h:by;
lambda = 0;
z = x'*y; %creates a N+2 X N+2 matrix to build a U matrix with starting BC. 

%iterative conditions
tol = 1e-6; %tolerance level
it = 10*N^2; %iterations scaling to N

%Dirichlet Boundary Conditions
uy1 = (bx - x).^2 .* cos(pi().*x./(bx)); %fb(x) top side BC
uy2 = x.*(bx -x).^2; %gb(x) bottom side BC
ux1 = (bx - ax).^2 * cos(pi()*ax/(bx)) + ((y - ay)/(by - ay))*(ax.*(bx-ax).^2 - (bx-ax).^2*cos(pi()*ax/(bx))); %Left Side Boundary Condition
%Neuman Boundary Condition
du_dx = 0;

%combing BC into a single matrix of size determined by N+2
U = zeros(size(z));
U(1,:) = uy1; 
U(N+2,:) = uy2; 
U(:,1) = ux1; 
U(2:N+1,N+2) = du_dx; 

%iterating for Given F equation
F = zeros();
for jj = 1:N+2
   for ii = 1:N+2
      F(ii,jj) = (sin(pi()*((x(ii)-ax)/(bx-ax))).*cos((pi()/2).*(2.*((y(jj)-ay)/(by-ay)) +1))); 
   end
end

%Approximation Using Gauss Seidel and SOR Methods
%when t = 1 Gauss Seidel is determined
%when 1>t>2 SOR is determined
w = 1.78; %Correction Factor for SOR
er = zeros();
error = .1 + tol;
iter = 0;
tic
while error > tol
    iter = iter + 1;
    Uold = U;
  for jj = 2:N+1
    for ii = 2:N+1
        U(ii,jj) = (1/(4-(lambda*h^2))*(U(ii-1,jj) + U(ii+1,jj) + U(ii,jj-1) + U(ii,jj+1) - (h^2*F(ii,jj))));
    end 
    %Solves Neuman Boundary
        U(N+2,jj) = (-2*U(N+1,jj) - U(N+2,jj-1) - U(N+2,jj+1) + h^2*(F(N+2,jj)))/(-4+(h^2*lambda));
   for ii = 2:N+1
          %Calculating error
     if U(ii,jj) ~= 0
        %Computing SOR when t is between 1 and 2. Otherwise it remains the same.
        U(ii,jj) = w*U(ii,jj)+(1-w)*Uold(ii,jj); 
        er(ii,jj) = abs((U(ii,jj) - Uold(ii,jj))/U(ii,jj))*100;
     end
   end
  end
         if iter > it
             break
         end
         error = max(max(er));
         if error <= tol
             fprintf('Tolerance achieved in %d iterations', iter);
             break
         end
end
toc
[xx,yy]=meshgrid(x',y);
    surf(xx,yy,U);
    title('Helmholtz Convergence')
    xlabel('x axis')
    ylabel('y axis')