p1=pwd;
p1=[p1(1:strfind(p1, 'My Dropbox')-1) 'My Dropbox\Data\Tomo\2011-06-21\2011-06-21_afterIntrefero-ABDFinterfero-xxxBox-oldRL-15s-01.txt'];
raw_counts = load(p1);
%Change the counting time according to the file name here!!!!!!!!!
counting_time=15;
coin_win = 190e-9;
MLHFlag=1;

projectors = zeros(16,4);
M = zeros(16,4,4);
global n;
n = zeros(16,1);
rho = zeros(4,4);



%projectors=proj_james;

projectors=proj_path1Interfero;

B = B_matrix(projectors);
B_inv = inv(B);

for mu=1:1:16
    M(mu,:,:) = M_matrix(mu, projectors, B, B_inv);
end

%if size(raw_counts,2)==5
    %coinc = raw_counts(:,5)'; % to disable/enable the filtering, set to column 4/5 here.
    %use the singles rate to calculate the accidentals and subtract.
    %coinc = raw_counts(:,4)-raw_counts(:,2).*raw_counts(:,3)*240e-9/counting_time*4;
    %coinc = coinc';
    %coinc = raw_counts(:,5)';
    %else
    
    %coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.95;
    %coinc = raw_counts(:,11) - raw_counts(:,7).*raw_counts(:,9)*coin_win/counting_time*4*0.92;
    %coinc = raw_counts(:,11)- raw_counts(:,7).*raw_counts(:,9)*coin_win/counting_time*4*1.015;
    %coinc = raw_counts(:,11) - raw_counts(:,7).*raw_counts(:,9)*coin_win/counting_time*4*0.6;    
    coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*0.92;
    
    %Correction factor 1.108, from accidental msmt. '2011-03-03_H_not
    %alighed. Blue 0.190v, delay 130ns, 15s, overcomplete basis - 01.txt
    
    %coinc = raw_counts(:,4) - 1/coin_win*(1 - exp(-raw_counts(:,2)/counting_time*4*coin_win)).*1/coin_win.*(1 - exp(-raw_counts(:,3)/counting_time*4*coin_win))...
    %                          * coin_win * counting_time / 4;

    %disp(coinc);
    
    ind1=find(coinc<0);
    coinc(ind1)=0.0001;
    coinc = coinc';
    %coinc = abs(coinc');
    %disp(coinc);

%end

if size(coinc,2)==36
    ind = [1; 2; 8; 7; 25; 26; 14; 13; 17; 15; 27; 3; 9; 12; 6; 30];
    n = coinc(ind);
else if size(coinc,2)==16
        n = coinc;
     end
end

%Data from Daniel's paper
% n = [34749; 324; 35805; 444; 16324; 17521; 13441; 16901; 17932; 32028; 15132; 17238; 13171; 17170; 16722; 33586];
%another fake data: theoretical
%n=[0; 1000; 0; 1000; 500; 500; 500; 500; 500; 1000; 500; 500; 500; 500; 500; 0];
%From the script: 'Counts_From_Rho.m'
%n=[0.00;500.00;0.00;500.00;250.00;250.00;250.00;250.00;250.00;0.00;250.00;250.00;250.00;250.00;250.00;500.00];

for nu=1:1:16
    rho = rho + reshape(M(nu,:,:),4,4)*n(nu);
end
rho = rho/sum(n(1:4));

% f1=figure;
% subplot(1,2,1);
% bar3(real(rho));
% set(gca,'XTickLabel',{'HH','HV','VH','VV'},'YTickLabel',{'HH','HV','VH','VV'})
% title('real part');
% 
% subplot(1,2,2);
% bar3(imag(rho));
% set(gca,'XTickLabel',{'HH','HV','VH','VV'},'YTickLabel',{'HH','HV','VH','VV'});
% title('imaginary part');
% 
% d1=eig(rho);
% purity1 = sum(diag(rho*rho));
% 
% rho_t = [0 0 0 0; 0 0.5 0.5*exp(-i*2.5) 0; 0 0.5*exp(i*2.5) 0.5 0;0 0 0 0];
% fidelity = sum(diag(rho*rho_t));
% 
% f2=figure; bar3(abs(rho))
% set(gca,'XTickLabel',{'HH','HV','VH','VV'},'YTickLabel',{'HH','HV','VH','VV'});
% title('Amplitude (Linear Inversion)');
% % figure; bar3(phase(rho))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Maximum likelyhood estimation
if MLHFlag
    if size(n,2)==36
        n= coinc; %use all the availible data for MLH reconstruction.
    end
    t=FindInitialT(rho);
    fhandle=@fun_MLH;
    
    [t,fval]=fminsearch(fhandle,t,optimset('MaxIter',1000*length(t),'MaxFunEvals',1000*length(t)));
    rho_mlh=fun_rho(t);
    
    d2=eig(rho_mlh);
    purity2 = sum(diag(rho_mlh*rho_mlh));
    
    figure; bar3(abs(rho_mlh))
    set(gca,'XTickLabel',{'HH','HV','VH','VV'},'YTickLabel',{'HH','HV','VH','VV'});
    title('Amplitude (MLH Estimation)');
    disp([,'Eig_Linear','   Eig_MLH']);
    disp(real([d1,d2]));
    disp(['Purity_Linear','  Purity_MLH']);
    disp(real([purity1,purity2]));
else
    disp('Eig_Linear');
    disp(real(d1));
    disp('Purity_Linear');
    disp(real(purity1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %find out the most probable state
% n=50;
% tht=linspace(0,2*pi,n); 
% ph=linspace(0,2*pi,n);
% %alp=linspace(0,2*pi,n);
% fidp=zeros(n,n);
% 
% temp=zeros(1,10);
% % alp=linspace(0,pi,10);
% alp=2*pi/3;
% % for indd=1:10
% for indi=1:n
%     for indj=1:n
% %        for indk=1:n
%             eh=[cos(tht(indi));sin(tht(indi))*exp(i*ph(indj))];
%             ev=[-sin(tht(indi))*exp(i*ph(indj));cos(tht(indi))];
% 
%             psip=tensor_product(eh,ev)+tensor_product(ev,eh)*exp(i*alp);
%             rhop=psip*psip'/2;
%             fidp(indi,indj) = sqrt(sum(diag(rho*rhop)));
% %        end
%     end
% end
% % temp(indd)=max(max(fidp));
% % end
% 
% figure;
% % plot(alp,temp)
% max(fidp)
% surf(tht,ph,abs(fidp));
% %colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
