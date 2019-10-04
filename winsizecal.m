function [sr,sc,er,ec]=winsizecal(rp,cp,N,dim,dim2)
difr=0;
difc=0;
if rp-N/2>0
    sr=rp-N/2;
else
    difr=N/2-rp;
    sr=1;
end
if ~(rp+N/2>dim)
    er=rp+N/2+difr;
else
    er=dim;
    difr=rp+N/2-dim;
    sr=sr-difr;
end
if cp-N/2>0
    sc=cp-N/2;
else
    sc=1;
    difc=N/2-cp;
end
if ~(cp+N/2>dim2)
    ec=cp+N/2+difc;
else
    ec=dim2;
    difc=cp+N/2-dim2;
    sc=sc-difc;
end
    
    
    
    
    