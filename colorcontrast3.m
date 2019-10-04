function [al1,au1,av1]= colorcontrast3(al_k,au_k,av_k,m,str)
% str=std(std(al));
[nn1,nn2]=size(al_k);
% M=max(max(al(:)));
al1(nn1,nn2)=0;
 au1(nn1,nn2)=0;
  av1(nn1,nn2)=0;
p=.01;
spo=0;
% m=.9*M;
no=4*str;
for kk1=1:nn1
    for kk2=1:nn2
        e=sqrt(al_k(kk1,kk2)^2+au_k(kk1,kk2)^2+av_k(kk1,kk2)^2);
%         e=abs(e);
         if e<no
            al1(kk1,kk2)=40*al_k(kk1,kk2);
             au1(kk1,kk2)=40*au_k(kk1,kk2);
              av1(kk1,kk2)=40*av_k(kk1,kk2);
        else
            if (no<e && e<(2*no))
%                 if (no<e <(2*no))
                s=((e-no)/no)*((m/no)^p)+(2*no-e)/no;
                s=2*abs(s);
                al1(kk1,kk2)=al_k(kk1,kk2)*s;
                 au1(kk1,kk2)=au_k(kk1,kk2)*s;
                  av1(kk1,kk2)=av_k(kk1,kk2)*s;
            else
%                 if  ((2*no)<e && e<m)
                      if  ((2*no)<e<m) 
                    s2=8*(m/e)^p;
                    al1(kk1,kk2)=s2*al_k(kk1,kk2);
                    au1(kk1,kk2)=s2*au_k(kk1,kk2);
                    av1(kk1,kk2)=s2*av_k(kk1,kk2);
                else
                    if  (m<e)
                         s3=1*(m/e)^spo;
                         al1(kk1,kk2)=s3*al_k(kk1,kk2);
                         au1(kk1,kk2)=s3*au_k(kk1,kk2);
                    av1(kk1,kk2)=s3*av_k(kk1,kk2);
                         
                    end
                end
            end
         end
    end
end