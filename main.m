function [imfi]= main(ju,refju,N,imn,s,mnc,iter,sigmaa,imna,agcoe,this_scan)
dim=size(imn,1);
dims2=size(imn,2);
imfi=zeros(size(imn));
w=zeros(size(imn));


hhh = waitbar(0,strcat('Please wait...','Iteration# ',num2str(iter),'....Scan# ',num2str(this_scan)));
steps = dim;
% for step = 1:steps
%     % computations take place here
%     waitbar(step / steps)
% end
% close(h) 


for rp=1:refju:size(imn,1)
    clc;
    iter
    waitbar(rp / steps)
    progress=rp/dim
    if (size(imn,1)-rp)<s
        rp=size(imn,1)-s+1;end
    for cp=1:refju:size(imn,2)
        if (size(imn,2)-cp)<s
            cp=size(imn,2)-s+1;end
        refp=imn(rp:rp+s-1,cp:cp+s-1);
        refpa=imna(rp:rp+s-1,cp:cp+s-1);
        refpv=reshape(refp,[s*s,1]);
        refpva=reshape(refpa,[s*s,1]);
        [sr,sc,er,ec]=winsizecal(rp,cp,N,dim,dims2);
        dimc=0;  
        cs=0;
        for opr=sr:ju:er-s+1
            for opc=sc:ju:ec-s+1
                cs=cs+1;
            end
        end
        Y1=zeros(s*s,cs);
        Y1a=zeros(s*s,cs);
        pos1=zeros(2,cs);
        for opr=sr:ju:er-s+1
            for opc=sc:ju:ec-s+1
                dimc=dimc+1;
                pa=imn(opr:opr+s-1,opc:opc+s-1);
                paa=imna(opr:opr+s-1,opc:opc+s-1);
                y=reshape(pa,[s*s,1]);
                ya=reshape(paa,[s*s,1]);
                Y1(:,dimc)=y;
                Y1a(:,dimc)=ya;
                pos1(:,dimc)=[opr;opc];
            end
        end
        Y=zeros(s*s,mnc+1);
        Ya=zeros(s*s,mnc+1);
        [Y,Ya,pos]=KNN(refpv,Y1,s,pos1,mnc,Y1a);
        dimc2=size(Y,2);
        Y(:,dimc2+1)=refpv;
        Ya(:,dimc2+1)=refpva;
        pos(:,dimc2+1)=[rp;cp];
        %denoising
        [Yd wc]=lowrcovsimp2(Y,Ya,sigmaa,agcoe);
        clear Y;
        for i=1:size(Yd,2)
            cor1=pos(1,i);
            cor2=pos(2,i);
            blk=reshape(Yd(:,i),[s,s]);
            imfi(cor1:cor1+s-1,cor2:cor2+s-1)=imfi(cor1:cor1+s-1,cor2:cor2+s-1)+blk*wc(i,1);
            w(cor1:cor1+s-1,cor2:cor2+s-1)=w(cor1:cor1+s-1,cor2:cor2+s-1)+wc(i,1);
        end
        
                   
    end
end
close(hhh)
%final measurements
imfi=imfi./w;