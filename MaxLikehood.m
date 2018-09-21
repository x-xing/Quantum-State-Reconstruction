clear; 
close all;

xtomo_linear;

t=FindInitialT(rho);

fhandle=@fun_MLH;

[t,fval]=fminsearch(fhandle,t,optimset('MaxIter',400*length(t),'MaxFunEvals',600*length(t)));

rho_mlh=fun_rho(t);

eig(rho_mlh)
purity = sum(diag(rho_mlh*rho_mlh))

figure; bar3(abs(rho_mlh))
set(gca,'XTickLabel',{'HH','HV','VH','VV'},'YTickLabel',{'HH','HV','VH','VV'});
title('Amplitude');