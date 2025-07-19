function Xnew = update_svd_thin(Xold, u_add, v_add)
%  update the svd decomposition for rank 1 updation matrix


Uold = Xold.U;
dold = Xold.d;
Vold = Xold.V;

if isempty(Uold) && isempty(Vold) && isempty(dold)
    u_norm = norm(u_add, 2);
    v_norm = norm(v_add, 2);
    Unew = u_add/u_norm;
    Vnew = v_add/v_norm;
    Dnew = u_norm*v_norm;
else
    [Unew, Dnew, Vnew] = svd_rank1_update_qr(Uold, diag(dold), Vold, u_add, v_add);
end

[Unew, Dnew, Vnew] = thinSVD(Unew, Dnew, Vnew, 1e-6);
Xnew = struct();
Xnew.U = Unew;
Xnew.V = Vnew;
Xnew.d = diag(Dnew);
