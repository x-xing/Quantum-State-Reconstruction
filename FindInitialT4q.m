function outT=FindInitialT4q(rho)
    rho = makephysical(rho);
    outT = chol(rho);
end
    