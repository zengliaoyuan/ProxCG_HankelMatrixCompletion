% This function takes in mi by nj matrix and output the Hankel adjoint: m
% by n(i+j-1) matrix

function Hkxad = Hkxad(H,m,n,i,j)

Hkxad = zeros(m,n*(i+j-1));
for k = 1:j
  if k <= i
    S = zeros(m,n*(i-k+1));
    for s = 1:(i-k+1)
      S(:,n*(s-1)+1:n*s) = H((s-1)*m+1:s*m,(k-1)*n+1:k*n);
    end
    Hkxad = Hkxad + [zeros(m,(k-1)*n) S H((i-k)*m+1:(i-k+1)*m,k*n+1:end) zeros(m,(k-1)*n)];
  end
end


