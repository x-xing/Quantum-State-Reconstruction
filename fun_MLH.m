function outL=fun_MLH(t,path,n,method)
%'path' parameter: 
%       0: normal tomography.
%       1: path tomo. 'uu'
%       2: path tomo. 'ul'
%       3: path tomo. 'lu'
%       4: path tomo. 'll'
%'method' paramter: 
%       0: using triangle T matrix, adopted from Daniel's paper.
%       2: using diagonalization to find T matrix.

Rho = fun_rho(t, method);

if size(n,2)==16
    projectors=proj_james;
    N=sum(n(1:4));
else if size(n,2)==36
        projectors=proj_full;
        N=sum([n(1:2) n(7:8)]);
    else
        error('Error in input count data. Neither 16 or 36!');
    end
end
switch path
    case 0
        ;
    case 1
        projectors = proj_2q_path('uu');
    case 2
        projectors = proj_2q_path('ul');
    case 3
        projectors = proj_2q_path('lu');
    case 4
        projectors = proj_2q_path('ll');
    otherwise
        error('Unknown ''path'' parameter');
end

dim_p = size(projectors, 1);
tmp=[];

for i=1:dim_p
    tmp = [tmp projectors(i,:)*Rho*(projectors(i,:)')];
end

outL=sum((N*tmp-n).^2/2./double(N*tmp));
