function [Unew, Dnew, Vnew]= svd_rank1_update_qr(Uold, Dold, Vold, u_add, v_add)

% % obtain svd of   Uold*Dold*Vold' + u_add*v_add' = Unew*Dnew*Vnew';
% % The method is from
% %"M. Brand. Fast low-rank modifications of the thin singular value decomposition, 2006."

[Q1, R1] = qr([Uold u_add], 0);
[Q2, R2] = qr([Vold v_add], 0);
z = zeros( size(R1, 2)-1, 1 );
K = R1*[ Dold z ; z' 1 ]*R2';

[tUp,tSp,tVp] = svd( K );


Dnew = tSp;

Unew = Q1 * tUp;
Vnew = Q2 * tVp;


