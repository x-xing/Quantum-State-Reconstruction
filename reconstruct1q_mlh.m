function out = reconstruct1q_mlh(rho_i,n_prob,proj)

rho=makephysical(rho_i);
t0 = chol(rho)';

options = optimset('MaxFunEvals',1e4,'MaxIter',1e4);
t=fminsearch(@fun1q_mlh,t0,options,n_prob,proj);
out=t*t';