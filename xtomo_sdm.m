%Comparison of two likelihood function.
%one is simple least square, the other is the correct likelihood which
%involves division of SDP variable. The latter was thought to be
%non-convex, but turns out to be fine with appropriate initialization, and
%using doulbe(SDP variable).
% Xingxing Xing, 2011-08-17

load('F:\xing\Dropbox\Data\Dylan''s code\counts.mat');
load('F:\xing\Dropbox\Data\Dylan''s code\M.mat');
load('F:\xing\Dropbox\Data\Dylan''s code\rho.mat');

np=zeros (1,36);
for ind=1:36
    np(ind)=trace(M(:,:,ind)*rho);
end
npe=[counts(:,1)' counts(:,2)' counts(:,3)' counts(:,4)'];
%plot(1:36,[np;npe]);

rhosdm=sdpvar(4,4,'hermitian','complex');
assign(rhosdm,rho);
val=0;
for ind=1:36
    n_sdm=trace(M(:,:,ind)*rhosdm);
    val=val+(n_sdm-npe(ind)).^2./double(n_sdm);
end
solvesdp([rhosdm>0,trace(rhosdm)==1],val);
rho_sdm = double(rhosdm);
rho_sdm=rho_sdm/trace(rho_sdm)
