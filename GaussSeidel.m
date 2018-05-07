function u = GaussSeidel(A,b,tol,it)
%A is the coefficient matrix for u
%b is the rhs vector
%tol is the tolerance, we will use 1e-6
%it is the iteration prescribed by the user, we will try up to 50 iterations
%u is the solution vector
[m,n] = size(A);
if m~=n, error('Matrix A must be square');end
   C = A;
for i = 1:n
    C(i,i) = 0;
    u(i) = 0;
end
u = u';
for i = 1:n
    C(i,1:n) = C(i,1:n)/A(i,i);
end
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
for i = 1:length(u)
    fprintf('\nu%d = %f\n',i,u(i));
end
end
    
