function [Y,Ya,pos]=KNN(ref,Y1,s,pos1,mc,Y1a)
dis1=zeros(1,size(Y1,2));

for i=1:size(Y1,2)
    dis1(i)=sum((Y1(:,i)-ref).^2)/s^2;
end

       [dis indx2]=sort(dis1);
        Y1=Y1(:,indx2);
        Y1a=Y1a(:,indx2);
        pos1=pos1(:,indx2);
        Y=Y1(:,1:mc);
        Ya=Y1a(:,1:mc);
        pos=pos1(:,1:mc);
        

        
    
    
    
    
