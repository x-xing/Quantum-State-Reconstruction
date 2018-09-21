out=zeros(2,180);
for i=1:1:180
    out(i)=QWP(i/180*pi)*PBS*[1;0];
end

plot(out)
    