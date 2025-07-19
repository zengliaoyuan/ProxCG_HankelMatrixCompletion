function x = Hkxad_fast(U,V)
% This computes the adjoint of the Hankel operator of UV', i.e.,  H^*(UV')
% The code is a part of the HSGD.m downloaded from
%  https://github.com/caesarcai/HSGD

n = size(U,1) + size(V, 1) -1;
r = size(V,2);
L = 2^nextpow2(n);
x = zeros(n,1);
for i = 1:r
    ui = U(:,i);
    vi = conj(V(:,i));
    ui = fft(ui,L);
    vi = fft(vi,L);
    ui = ui.*vi;
    ui = ifft(ui);
    ui = ui(1:n);
    x = x+ui;
end