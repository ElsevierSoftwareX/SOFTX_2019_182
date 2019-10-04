function resfast = fastdirec( potential,QQ,ii )

potential=potential;
[xx yy]=find(QQ~=0);
start_point=[xx';yy'];
clear options;
options.nb_iter_max = Inf;
options.end_points = [];
[D,S,numiter] = perform_fastmarching(potential, start_point, options);
[xx1 yy1]=find(D~=inf & D~=0);
test=zeros(1,length(xx1));
for i=1:length(test)
    test(1,i)=D(xx1(i),yy1(i));end
me_test=mean(test);std_test=std(test);
[xx2 yy2]=find(test<me_test);
test1=zeros(1,length(xx2));
for i=1:length(xx2)
    test1(1,i)=test(xx2(i),yy2(i));end
me_test1=mean(test1);std_test1=std(test1);
[xx3,yy3]=find(D==0 | D<abs(me_test1-std_test1));
resfast=zeros(size(D));
for i=1:length(xx3)
    resfast(xx3(i),yy3(i))=ii;end
