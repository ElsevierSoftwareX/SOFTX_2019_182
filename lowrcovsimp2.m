function [Yd we]=lowrcovsimp2(Y,Ya,sigmaa,agcoe)

  Yd=zeros(size(Y));
    las=1;
    ep=.1;
    ym=mean(Y,2);
    for iter=1:1
        dimy=size(Y,2);
        Sc=zeros(size(Y,1),size(Y,1));
    for i=1:dimy
        y=Y(:,i);
        y=y-ym;
        covy=(y)*(y)';
        Sc=covy+Sc;
    end
    covY=Sc./(dimy-1);

    [U S]=eig(covY);
    Si=zeros(size(S));
    m=size(S,1);
    for i=1:m
            S(i,i)=S(i,i)-0;
            Si(i,i)=1/(S(i,i)+ep);
    end
    Ym=zeros(size(Y));
    S=S+ep*eye(size(S));
    for i=1:size(Y,2)
        Ym(:,i)=ym;
    end
    X=U'*Ya+1*(las*(sigmaa^2)*Si)*U'*Ym;
    for i=1:size(S,1)
    X(i,:)=(1/(1+las*(sigmaa^2)*Si(i,i)))*X(i,:);
    end  
    Yd=U*X;
    dets=1;
    for i=1:m
        dets=dets*S(i,i);
    end
    we=ones(size(Y,2),1);
    for i=1:size(Y,2)
        l=((Yd(:,i)-Ya(:,i))'*(Yd(:,i)-Ya(:,i))+las*sigmaa^2*(Yd(:,i)-ym)'*(U*Si*U')*(Yd(:,i)-ym)+las*sigmaa^2*log(dets))*agcoe;
        we(i,1)=exp(-l);
    end
        


    end