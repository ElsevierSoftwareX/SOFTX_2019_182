function D= edgenhance(D)
str=std(std(D));
[n1,n2]=size(D);
M2=max(max(D(:)));
N2=min(min(D(:)));
for k1=1:n1
    for k2=1:n2
        if D(k1,k2)>0
%             D(k1,k2)=D(k1,k2)*2;N2*str+(M2-N2)
D(k1,k2)=0;
        else
%              D(k1,k2)= D(k1,k2)*7;
D(k1,k2)=0;
        end
    end
end
