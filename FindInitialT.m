function outT=FindInitialT(rho,method)
%method =0: Directly adopted from Daniel's paper. The analytical form is used.
%method =1: Using diagonalization to get t.

if ~method
    outT=zeros(1,16);
    outT(1)   =   real(sqrt(det(rho)/minor(rho,1,1)));
    temp      =   minor(rho,1,2)/sqrt(minor(rho,1,1)*minor(rho,1:2,1:2));
    outT(2)   =   real(temp);
    outT(3)   =   imag(temp);
    outT(4)   =   real(sqrt(minor(rho,1,1)/minor(rho,1:2,1:2)));
    temp      =   minor(rho,1:2,2:3)/sqrt(rho(4,4))/sqrt(minor(rho,1:2,1:2));
    outT(5)   =   real(temp);
    outT(6)   =   imag(temp);
    temp      =   minor(rho,1:2,[1 3])/sqrt(rho(4,4))/sqrt(minor(rho,1:2,1:2));
    outT(7)   =   real(temp);
    outT(8)   =   imag(temp);
    outT(9)   =   real(sqrt(minor(rho,1:2,1:2)/rho(4,4)));
    temp      =   rho(4,1)/sqrt(rho(4,4));
    outT(10)  =   real(temp);
    outT(11)  =   imag(temp);
    temp      =   rho(4,2)/sqrt(rho(4,4));
    outT(12)  =   real(temp);
    outT(13)  =   imag(temp);
    temp      =   rho(4,3)/sqrt(rho(4,4));
    outT(14)  =   real(temp);
    outT(15)  =   imag(temp);
    outT(16)  =   real(sqrt(rho(4,4)));
else if method == 1
        rho=makephysical(rho);
        outT = chol(rho)';
    else
        error('unknown method parameter in FindInitialT.');
    end
end


