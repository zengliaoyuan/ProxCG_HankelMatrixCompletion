function [proj_sy, flag] = project_sigVal_to_simplex(s, sigma)
% The s should be nonnegative and in a descent order. 
% The output of the Matlab svd is nonnegative and in descent order.
    flag = false;
    proj_sy = s;
    dim_s = size(s);
    s0= 0;
    for i = 1: dim_s(1)
        s0 = s0+ s(i);
        if s0 - sigma >=i*s(i)
            s0 = s0 - s(i);
             number_s = i-1;
            break;
        else
            number_s = i;
        end

    end
    if s0 > sigma
        proj_sy = max(s- (s0 - sigma)/number_s,0);
        flag = true;
    end
end