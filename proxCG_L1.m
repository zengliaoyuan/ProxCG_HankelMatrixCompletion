function [x, y, fval, iter, cpu_time, history]  = proxCG_L1(w, obs, Omega, sigma, m, n, opts)
%  This solves
% min ||w.*(x+s-obs)||_1
% s.t. ||H(x)||_* <= sigma,
% x is a m+n-1 by 1 complex vector,
% H(x) is a m by n Hankel matrix formed from x,
% w is the weight matrix, and obs come from measurements.

% % We always initial y0 = H(x0) and we do not form y0 explicitly

if isfield(opts, 'maxiter')
    maxiter = opts.maxiter;
else
    maxiter = inf;
end

if isfield(opts, 'scale_beta')
    scale_beta = opts.scale_beta;
else
    scale_beta = 1;
end

cputime_now = 0;
tstart = tic;

if isfield(opts, 'xinit')
    x = opts.xinit;
else
    x = zeros(m+n-1, 1);
    x(Omega) = obs;
end

% ss = inf;
lambda_H = min(m, n);     % || H^*H ||
w_Omg = w(Omega);
sigmatilde = sigma + norm(obs) + 1;
H_0 = 1e-6;

fval_rec = norm(w_Omg.*(x(Omega) -obs),1);
ranky_rec = NaN;    % We do not form y0 explicitly and we have not the infomation of  rank(y0)
iter = 0;
%% main loop
while 1
    
    % update x
    xold = x;
    
    beta_t = sqrt(iter+1)* scale_beta;
    if iter == 0
        x = xold;        % x1 = x0;
    else
        HadjHx = w.*x;
        Hadjy = Hkxad_fast(y.d'.*y.U, y.V);
        Hadj_yHx = HadjHx - Hadjy ;
        
        coefficient = H_0 + beta_t * lambda_H;
        x = x - beta_t/coefficient * Hadj_yHx;
        prox_const = x(Omega) - obs;
        z = x;
        z(Omega) =  sign(prox_const).*max(abs(prox_const) - w_Omg / coefficient, 0);            % z = x - obs;
        
        nm_z = norm(z);
        if nm_z  <= sigmatilde
            x(Omega) = obs + z(Omega);
        else
            z_scale = z*(sigmatilde/nm_z);
            x = z_scale;
            x(Omega) = obs + x(Omega);
        end
    end
    
    fval = norm(w_Omg.*(x(Omega) -obs),1);
    fval_rec = [fval_rec; fval];
    
    %Update y
    alpha = 2/(iter+2);
    if iter == 0
        ss = 0;    
        y.d = 0;        % note that alpha0 = 1 and y0 - H(x1) = H(x0) - H(x1)=0. So we have y1 = 0;
        y.U = zeros(m, 1);
        y.V = zeros(n, 1);
    else
        [u, ss, v] = svds(@(w, tflag)yHx_hdl(w, tflag, y.U, y.d, y.V, x),  [m, n], 1);
        ybase = y;
        ybase.d = (1- alpha)*y.d;
        y = update_svd_thin(ybase, sigma*u, -alpha*v);
    end
    ranky_rec = [ranky_rec; sum(y.d > 1e-6)];
    
    % % display
    if mod(iter, opts.verbose_freq) == 0
        fprintf('iter No. %7d Time %7.1f: fval = %6.4e, num_rank(y) = %d, rel_nuc_feas = %2.1e, rel. ||Y-Hx_+|| = %2.1e \n',...
            iter+1, toc(tstart), fval, sum(y.d/max(max(y.d), 1) > 1e-5), sum(y.d)/sigma-1, ss/max(max(y.d), 1) )
        if (sum(y.d)-sigma) > 1e-12
            fprintf('y infeasible: ||y||_* - sigma = %6.2e\n', sum(y.d)-sigma)
        end
    end
    
    % %     termination
%     if norm(x-xold)/max(norm(xold), 1) <= 1e-12  && iter > 0     % x1 = x0
%         fprintf('Termination: the succesive changes is small. \n');
%         iter = iter + 1;    
%         break;
%     end
    
    if iter > maxiter-1   % note that iter starts from 0
        fprintf('Termination: maximal iteration number used.\n');
        iter = iter + 1;    
        break;
    end
    
    %% cputime
    iter = iter+1;
    cputime_now = toc(tstart);
    
end

%%
cpu_time = cputime_now;
history.fval_rec = fval_rec;
history.rank_rec = ranky_rec;
