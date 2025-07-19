function [x, y, fval, iter, cpu_time, history]  = completeADM(w, obs, Omega, sigma, m, n, opts)
%  This solves
% min ||w.*(x-obs)||_1
% s.t. ||H(x)||_* <= sigma,
% x is a m+n-1 by 1 complex vector,
% H(x) is a m by n Hankel matrix formed from x,
% w is the weight matrix, and obs come from measurements.

if isfield(opts, 'maxiter')
    maxiter = opts.maxiter;
else
    maxiter = inf;
end


fval_rec = [];
% lc_res_rec = zeros(1,1);

cputime_now = 0;
tstart = tic;

w_Omg = w(Omega);


if isfield(opts, 'y_ADMMinit')
    y = opts.y_ADMMinit;   
else
    y = zeros(m,n);
end

if isfield(opts, 'multiplier')
    Lambda = opts.multiplier;
    norm_lam = norm(Lambda, 2);
else
    Lambda  = zeros(m,n);
    norm_lam = 0;
end

iter = 1;
beta = 1;

Hadj_Lam = Hkxad(Lambda, 1,1,m,n).';
Hadjy = Hkxad(y,1,1,m,n).';
while 1  
        
    % update x
    winv = 1./ w;
    x_circ = -1/beta*winv.*Hadj_Lam  + winv.*Hadjy;
    
    vec_prox = -obs + x_circ(Omega);
    x_circ(Omega) = obs + sign(vec_prox).*max(abs(vec_prox) - 1/ beta, 0);
    x = x_circ;

    fval = norm(w_Omg.*(x(Omega) -obs),1);
    fval_rec = [fval_rec; fval];
    
    %update y.
    H_x = Hkx(x.',1,1,m,n);
    
    yhat = H_x+ 1/ beta*Lambda;
    
    [U_y, s_y, V_y] = svd(yhat,'econ');
    
    % Project the vector of the singular values onto a simplex.
    sy_vec = diag(s_y);
    [proj_sy_vec, flag] = project_sigVal_to_simplex(sy_vec, sigma);
    
    if flag
        y = proj_sy_vec'.*U_y*V_y';
    end 
    Hadjy = Hkxad(y,1,1,m,n).';

    
    %update multiplier
    y_minus_Hx = y - H_x;
    Lambda = Lambda - beta*y_minus_Hx;

    Hadj_Lam = Hkxad(Lambda, 1,1,m,n).';
    norm_lam = beta*max(sy_vec - proj_sy_vec);

    %compute the dual function value and dual feasibility
    Hadj_lam_omg = Hadj_Lam(Omega);
    dual_val = real(obs'*Hadj_lam_omg) - sigma*norm_lam;
    Hadj_lam_temp = Hadj_Lam;
    Hadj_lam_temp(Omega) = 0;
    dual_feas = norm(Hadj_lam_temp,1) + sum(max(abs(Hadj_Lam(Omega))- w_Omg,0));
    
    
    % % display
    if mod(iter, opts.verbose_freq) == 1
        fprintf('iteration %d  Time %7.1f： function_val = %6.4e, dual_val=%6.4e, rel_gap = %6.4e, rel_dual_feas = %2.1e\n',iter,...
            toc(tstart), fval, dual_val, abs(fval - dual_val)/fval, dual_feas/ norm_lam)
    end
    
    % termination 
        
    %  terminate if dual gap and the violation of dual feasibility is small
    if abs(fval - dual_val)/max(fval, 1) < 1e-1 && dual_feas/max(norm_lam, 1) < 5e-2
        fprintf('Termination: primal dual gap is small enough.\n');
        break;
    end
    
    % terminate if maximal iterations used
    if iter > maxiter
        fprintf('Termination: maximal iteration number used.\n');
        break;
    end
    
    %% cputime
    iter = iter+1;
    cputime_now = toc(tstart);
    
end
cpu_time = cputime_now;
history.fval_rec = fval_rec;