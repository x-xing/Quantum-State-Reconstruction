
 x=linspace(-30,30,100);
 y1=1./pi./(7^2+x.^2);
 o=zeros(1,100);
 
 for i=1:100
     delta=i*0.2;
     y2=1./pi./(7^2+(x-delta).^2);
 %the overlap
     o(i)=sum(y1.*y2)/sum(y1.*y1);
 end
 
 plot((1:100)*0.2,o)