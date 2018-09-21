function rhoMatrix=fun_rho(t,method)
%method = 0: from Daniel's paper. The triangle form guarantees hermition of rho.
%method = 1: using Cholesky factorization.
if ~method
tMatrix = [t(1)          0             0             0;...
           t(2)+i*t(3)   t(4)          0             0;...
           t(5)+i*t(6)   t(7)+i*t(8)   t(9)          0;...
           t(10)+i*t(11) t(12)+i*t(13) t(14)+i*t(15) t(16)];
rhoMatrix=tMatrix*tMatrix';

else if method == 1
        rhoMatrix = t*t';
    else
        error('unknown method in fun_rho');
    end
end
rhoMatrix=rhoMatrix/trace(rhoMatrix);