function [u,err] = GaussSeidel(A,b,tol,it)
%A is the coefficient matrix for u
%b is the rhs vector
%tol is the tolerance, we will use 1e-6
%it is the iteration prescribed by the user, we will try up to 50 iterations
%u is the solution vector
[m,n] = size(A);
if m~=n, error('Matrix A must be square');end
   C = A;
   u = zeros();
   if prod(diag(A))==0
disp('gauss-seidel iterative method is not applicable');
return
   end
   %Convergence Check
E=tril(A,-1);
F=triu(A,1);
D=A-E-F;
Bj=-inv(D)*(E+F);
rho=max(abs(eig(Bj))); %spectral radius
if rho>=1
disp('gauss-seidel method do not converge');
return
end
for i = 1:n
    C(i,i) = 0;
    u(i) = 0;
end
u = u';
for i = 1:n
    C(i,1:n) = C(i,1:n)/A(i,i);
end
d = zeros();
for i = 1:n
    d(i) = b(i)/A(i,i);
end
iter = 0;
while (1)
    uold = u;
    for i = 1:n
        u(i) = d(i)-C(i,:)*u;
        if u(i) ~= 0 
            err(i) = abs((u(i) - uold(i))/u(i)) * 100; %Calculates Normalized Error
        end
    end
    iter = iter + 1;
    if max(err)<= tol || iter >= it, break, end  %checks error against tolerance, or if max 'iter' is reached
end
end
    
