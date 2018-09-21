f1=figure;
subplot(1,2,1);
bar3(real(rho_mlh));
set(gca,'XTickLabel',{'HH','HV','VH','VV'},'YTickLabel',{'HH','HV','VH','VV'})
title('real part');

subplot(1,2,2);
bar3(imag(rho_mlh));
set(gca,'XTickLabel',{'HH','HV','VH','VV'},'YTickLabel',{'HH','HV','VH','VV'});
title('imaginary part');