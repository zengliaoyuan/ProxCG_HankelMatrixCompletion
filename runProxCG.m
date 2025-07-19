clear all
addpath('./PROPACK')

rng(1237)

n_val = [2^10; 2^12; 2^14; 2^16];
% n_val = [2^10; 2^12; 2^14];
r = 7;

reps = 10;


fprintf(' Start \n')
a = clock;
fname = ['run-table'  '-'  date  '-'  int2str(a(4))  '-'  int2str(a(5))  '.txt'];
fid = fopen(fname, 'w');
fprintf(fid, '%8s & %8s & %8s  & %8s  & %8s & %8s & %10s & %10s & %10s \n', ...
    'method',   'size', 'err',  'obj', 'iter', 'cpu' ,  'rel. $\sigma_r$', 'res. $\sigma_{r+1}$', 'rel.$\|\mathcal{H}(x_{out})\|$');

for nidx = 1: length(n_val)
    n = n_val(nidx);
    
    data_pCG = [];
    data_ADM = [];
    
    sample_rate = 0.4; % sample rate
    
    for repeat = 1:reps
        m = round(sample_rate*n);
        [Omega,ox,~] = generate_signal(n,r,m,'true','false');
        
        ox_orig = ox; % save the original signal
        
        noise_std = 1e-2;
        % Laplace noise;
        noise_real = randn(size(ox)).*randn(size(ox))- randn(size(ox)).*randn(size(ox));
        noise_img = randn(size(ox)).*randn(size(ox))- randn(size(ox)).*randn(size(ox));
        laplace_noise = noise_std*(noise_real+ noise_img*1i);
        
        ox = ox + laplace_noise;
        obs = ox(Omega);
        
        if mod(n,2)
            p = (n+1)/2;  %not sample rate here
            w = [1:p p-1:-1:1].';
        else
            p = n/2;
            w = [1:p p p-1:-1:1].';
        end
        q = n+1-p;
        
        %         %use lansvd
        %         Yforward = @(y) fhmvmultiply_1D(ox_orig,y);
        %         Ytranspose = @(y) fhmvmultiply_1D(conj(ox_orig),y);
        %         opts_svd = []; opts_svd.eta = 1e-16;
        %
        %         [~,ss_ox,~] = lansvd(Yforward,Ytranspose,p,q,r,'L',opts_svd);
        [~, ss_ox, ~] = svds(@(w, tflag)Hx_hdl(w, tflag, ox_orig),  [p, q], r);
        sigma = 0.97*sum(diag(ss_ox));
        % Assumed a good knowledge of the original signal. Likely the most critical parameter.
        
        fprintf('n = %d: The %d-th instances.\n\n', n, repeat)
        
        %%  ProxCG
        fprintf('Start of the ProxCG method. \n')
        % initial y and x, assuming explicit knowledge of r
        x_init = zeros(n,1);
        x_init(Omega) = obs;           % x_init = x_obs;
        opts.xinit = x_init;         % We always initial y0 as H(x_init) in proxFW_L1. So we do not set y0 here.
        
        opts.maxiter = 50000;
        opts.verbose_freq = 1000;
        opts.scale_beta = 0.3;
        [x_pCG, y_pCG, fval_pCG, iter_pCG, cpu_time_pCG, history_pCG]  = proxCG_L1(w, obs, Omega, sigma, p, q, opts);
        rec_err_pCG = norm(ox_orig - x_pCG)/norm(ox_orig);
        
        
        [~,s_Hx_pCG,~] = svds(@(w, tflag)Hx_hdl(w, tflag, x_pCG),  [p, q], r+1);
        s_Hx_pCG = diag(s_Hx_pCG);
        rel_sv_pCG = [s_Hx_pCG(r)/ s_Hx_pCG(1), s_Hx_pCG(r+1)/s_Hx_pCG(1)];
        rel_feas_pCG = sum(s_Hx_pCG)/sigma - 1;         % relative feasibility
        
        % save the results
        data_pCG = [data_pCG;  rec_err_pCG, fval_pCG, iter_pCG, cpu_time_pCG, rel_sv_pCG(1), rel_sv_pCG(2), rel_feas_pCG ];
        fprintf('proxCG end: rec_err = %6.4f, fval = %6.2e, iter = %d, time = %8.2f, rel. s_r = %6.2e, s_{r+1}=%6.2e, feas = %6.2e \n\n', ...
            rec_err_pCG, fval_pCG, iter_pCG, cpu_time_pCG, rel_sv_pCG(1), rel_sv_pCG(2), rel_feas_pCG);
        
        
        %% ADMM
        if n ~= 2^16
            fprintf('Start of ADMM. \n')
            % parameters for ADMM
            opts2.maxiter = 10000;
            opts2.verbose_freq = 200;
            opts2.y_ADMMinit = Hkx(x_init.', 1, 1, p, q);
            opts2.multiplier = zeros(p, q);
            
            [x_ADMM, y_ADMM, fval_ADM, iter_ADM, cpu_time_ADM, history_ADM]  = completeADM(w, obs, Omega, sigma, p, q, opts2);
            rec_err_ADMM = norm(ox_orig-x_ADMM)/norm(ox_orig);
            
            [~,s_Hx_ADM,~] = svds(@(w, tflag)Hx_hdl(w, tflag, x_ADMM),  [p, q], r+1);
            s_Hx_ADM = diag(s_Hx_ADM);
            rel_sv_ADM = [s_Hx_ADM(r)/ s_Hx_ADM(1), s_Hx_ADM(r+1)/s_Hx_ADM(1)];
            rel_feas_ADM = sum(s_Hx_ADM)/sigma - 1;         % relative feasibility
            
            
            % Save the results
            data_ADM = [data_ADM;  rec_err_ADMM, fval_ADM, iter_ADM, cpu_time_ADM, rel_sv_ADM(1), rel_sv_ADM(2), rel_feas_ADM ];
            fprintf('ADMM end: rec_err = %6.4f, fval = %6.2e, iter = %d, time = %8.2f, rel. s_r = %6.2e, s_{r+1}=%6.2e, feas = %6.2e \n\n', ...
                rec_err_ADMM, fval_ADM, iter_ADM, cpu_time_ADM, rel_sv_ADM(1), rel_sv_ADM(2), rel_feas_ADM);
        end
        
    end
    formatStrings = {'instance no.'; 'RecoveryError '; 'fval ';  'iter '; 'CPU time '; 'sigma_r_rel';  'sigma_{r+1}_rel:';  'rel_nuc_feas:';};
    Table_pCG = table(formatStrings, [1:reps; data_pCG']);
    if n ~=  2^16
        Table_ADM = table(formatStrings, [1:reps; data_ADM']);
        data_n = {Table_pCG, Table_ADM};
    else
        data_n = {Table_pCG};
    end
    dataname = ['Data_' 'n-' num2str(n) '_date' date '-' int2str(a(4)) '-' int2str(a(5)) '.mat'];
    save(dataname, 'n', 'data_n');    % save the result of every tested instance with dimension n
    
    fprintf(fid,  '%8s  &  %8d & %8.4f  & %8.2e  & %8.2f & %8.2f & %6.2e & %6.2e & %6.2e \n', ...
        'proxCG', n, mean(data_pCG));
    if n~=2^16
        fprintf(fid, '%8s  &  %8d & %8.4f  & %8.2e  & %8.2f & %8.2f & %6.2e & %6.2e & %6.2e \n', ...
            'ADMM', n, mean(data_ADM));
    end
end

fclose(fid);


function z = Hx_hdl(w, tflag, x)
% the function handle of H(x)w and conj(Hx)w

if strcmp(tflag,'notransp')
    Hxw = fhmvmultiply_1D(x, w);
    z = Hxw;
else
    Hxtw = fhmvmultiply_1D(conj(x),w);
    z = Hxtw;
end
end

