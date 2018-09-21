p1=pwd;
p1=[p1(1:strfind(p1, 'Dropbox')-1) 'Dropbox\Data\Tomo\2009-08-20\2009-08-20_tomo_10s, locked-Rb line1 singles10000-11.txt'];
raw_counts = load(p1);

%raw_counts = load('C:\\xingxing\\My Dropbox\\Data\\Tomo\\2011-08-08\\InterferoCorrected-maxOnDD-Mclose0far0-15s-01.txt');
%Change the counting time according to the file name here!!!!!!!!!
counting_time=10;
coin_win = 190e-9;
M = zeros(16,4,4);
n = zeros(16,1);
rho = zeros(4,4);
MLHFlag=1; %0 - no MLH; 1- MLH with fminsearh from MATLAB; 2-MLH with Sedumi.
method=0;  %0 - Analytic form of T matrix from Daniel's paper; 1 - from Cholesky factorization.
ftype = 0; %tomography mode: 0, normal; 1-4: path-tomography, 'uu','ul','lu','ll'.
sigmax=[0 1;1 0];


projectors = proj_james;

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
    %coinc = raw_counts(:,11) - raw_counts(:,7).*raw_counts(:,9)*coin_win/counting_time*4;
    %coinc = raw_counts(:,11)- raw_counts(:,7).*raw_counts(:,9)*coin_win/counting_time*4*1.015;
    %coinc = raw_counts(:,11) - raw_counts(:,7).*raw_counts(:,9)*coin_win/counting_time*4*0.6;    
    coinc = raw_counts(:,4) - raw_counts(:,2).*raw_counts(:,3)*coin_win/counting_time*4*1.04;
    
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

%%%% Diagonosis for 4 qubit tomography %%%%ss
% load('N_256Counts');
% n=n(16:16:256);

%Data from Daniel's paper
%n = [34749; 324; 35805; 444; 16324; 17521; 13441; 16901; 17932; 32028; 15132; 17238; 13171; 17170; 16722; 33586];
%another fake data: theoretical
%n=[0; 1000; 0; 1000; 500; 500; 500; 500; 500; 1000; 500; 500; 500; 500; 500; 0];
%From the script: 'Counts_From_Rho.m'
%n=[0.00;500.00;0.00;500.00;250.00;250.00;250.00;250.00;250.00;0.00;250.00;250.00;250.00;250.00;250.00;500.00];
%Dylan's counts probability to check Sedumi
%n=[];
drho = zeros(4,4);
for nu=1:1:16
    rho = rho + reshape(M(nu,:,:),4,4)*n(nu);
    drho = sqrt(drho.^2+(reshape(M(nu,:,:),4,4)*sqrt(n(nu))).^2);
end
drho= rho/sum(n(1:4)).*sqrt((drho./rho).^2+1/sum(n(1:4)));
rho = rho/sum(n(1:4));
switch ftype
    case 0
    case 1
        rho = tensor_product(sigmax,sigmax)*rho*tensor_product(sigmax,sigmax)';
        disp('Sigmax operation applied.');
%     case 2
%         projectors=proj_path('ul');
%     case 3
%         projectors=proj_path('lu');
%     case 4
%         projectors=proj_path('ll');
    otherwise
        warning('unknown tomography mode: ftype.');
end

figure(1);
subplot(1,2,1);
bar3(real(rho));
set(gca,'XTickLabel',{'HH','HV','VH','VV'},'YTickLabel',{'HH','HV','VH','VV'})
title('real part');

subplot(1,2,2);
bar3(imag(rho));
set(gca,'XTickLabel',{'HH','HV','VH','VV'},'YTickLabel',{'HH','HV','VH','VV'});
title('imaginary part');

d1=eig(rho);
purity1 = sum(diag(rho*rho));

rho_t = [0 0 0 0; 0 0.5 0.5*exp(-i*2.5) 0; 0 0.5*exp(i*2.5) 0.5 0;0 0 0 0];
fidelity(rho,rho_t);

figure(2); bar3(abs(rho))
set(gca,'XTickLabel',{'HH','HV','VH','VV'},'YTickLabel',{'HH','HV','VH','VV'});
title('Amplitude (Linear Inversion)');
% figure; bar3(phase(rho))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Maximum likelyhood estimation
if MLHFlag ==1
    if size(n,2)==36
        n= coinc; %use all the availible data for MLH reconstruction.
    end
    t0=FindInitialT(rho,method);
    options = optimset('MaxFunEvals',2e4,'MaxIter',1e4);
    [t,fval]=fminsearch(@fun_MLH,t0,options,ftype,n,method);
    rho_mlh=fun_rho(t,method);
    
    d2=eig(rho_mlh);
    purity2 = sum(diag(rho_mlh*rho_mlh));
    
    figure(3); bar3(abs(rho_mlh))
    set(gca,'XTickLabel',{'HH','HV','VH','VV'},'YTickLabel',{'HH','HV','VH','VV'});
    title('Amplitude (MLH Estimation)');
    disp(['Eig_Linear','   Eig_MLH']);
    disp(real([d1,d2]));
    disp(['Purity_Linear','  Purity_MLH']);
    disp(real([purity1,purity2]));
else if MLHFlag ==0
        rho
        disp('Eig_Linear');
        disp(real(d1));
        disp('Purity_Linear');
        disp(real(purity1));
    else if MLHFlag ==2
            %Sedumi solver
            if size(n,2)==36
                n= coinc; %use all the availible data for MLH reconstruction.
                projectors=proj_full;
                norm1 = sum(n(1:2))+sum(n(7:8));
            else if size(n,2)==16
                norm1 = sum(n(1:4));
                end
            end
            
            
            rhosdm=sdpvar(4,4,'hermitian','complex');
            assign(rhosdm,rho);
            val=0;
            for ind=1:16
                n_sdm=projectors(ind,:)*rhosdm*(projectors(ind,:)');
                val=val+(n_sdm-n(ind)/norm1).^2;%./double(n_sdm);
            end
            solvesdp([rhosdm>0,trace(rhosdm)==1],val);
            rho_sdm = double(rhosdm);
            rho_sdm=rho_sdm/trace(rho_sdm);
            d2=eig(rho_sdm);
            purity2 = trace(rho_sdm^2);
            figure(3); bar3(abs(rho_sdm))
            set(gca,'XTickLabel',{'HH','HV','VH','VV'},'YTickLabel',{'HH','HV','VH','VV'});
            title('Amplitude (MLH Estimation)');
            disp(['Eig_Linear','   Eig_SDM_MLH']);
            disp(real([d1,d2]));
            disp(['Purity_Linear','  Purity_SDM_MLH']);
            disp(real([purity1,purity2]));
        end
    end
end


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
%             fidp(indi,indj) = (sum(diag(rho*rhop)))^0.5;
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
 
