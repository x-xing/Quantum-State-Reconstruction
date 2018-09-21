a=0.5; %superposition parameter D=sqrt(1-a^2)H+aV.
n_prob = [1 0 0.75 .25 .75 .25];
proj = proj1q(a);
n_prob_c = n_prob([1 2 3 5]);
proj_c = proj([1 2 3 5],:);

rho_i = reconstruct1q_li(n_prob_c, proj_c);
rho_m = reconstruct1q_mlh(rho_i,n_prob,proj);

% n_prob_t = zeros(size(n_prob));
% for ind = 1:size(proj,1)
%     n_prob_t(ind)=trace(proj(ind,:)*rho_m*proj(ind,:)');
% end
% n_prob
% n_prob_t

% rho_t = [0.4 0.1+0.2i;0.1-0.2i 0.6];
% n_prob_t = zeros(size(n_prob));
% for ind = 1:size(proj,1)
%     n_prob_t(ind)=trace(proj(ind,:)*rho_t*proj(ind,:)');
% end
% n_prob_c = n_prob_t([1 2 3 5]);
% proj_c = proj([1 2 3 5],:);
% rho_i = reconstruct1q_li(n_prob_c, proj_c);
% rho_m = reconstruct1q_mlh(rho_i,n_prob_t,proj);
% rho_t
% rho_m