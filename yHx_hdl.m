function z = yHx_hdl(w, tflag, U, d, V, x)

% the function handle of (y-Hx)w  and conj(y-Hx)w

if strcmp(tflag,'notransp')
    Hxw = fhmvmultiply_1D(x, w);
    z = (d'.*U)*(V'*w) - Hxw;    % y = (UdV')*w - Hx*w;
else
    Hxtw = fhmvmultiply_1D(conj(x),w);
    z = (d'.*V)*(U'*w) - Hxtw;
end



