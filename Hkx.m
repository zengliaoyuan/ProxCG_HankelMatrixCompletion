% This function takes in an m by n(i+j-1) matrix X and output a block Hankel
% matrix of size mi by nj.

function Hkx = Hkx(X,m,n,i,j)

Hkx = zeros(m*i,n*j);

for k = 1:i
  Hkx((k-1)*m+1:k*m,1:n) = X(:,(k-1)*n+1:k*n);
end

for k = 2:j
  Hkx(1:(i-1)*m,(k-1)*n+1:k*n) = Hkx(m+1:i*m,(k-2)*n+1:(k-1)*n);
  Hkx((i-1)*m+1:i*m,(k-1)*n+1:k*n) = X(:,(i+k-2)*n+1:(i+k-1)*n);
end

