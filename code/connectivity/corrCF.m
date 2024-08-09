function [ro]=corrCF(q1,q2)

qq1=q1-mean(mean(q1));
qq2=q2-mean(mean(q2));
n1=sqrt(mean(mean(qq1.*qq1)));
n2=sqrt(mean(mean(qq2.*qq2)));

ro=mean(mean(qq1.*qq2))/(n1*n2);

end