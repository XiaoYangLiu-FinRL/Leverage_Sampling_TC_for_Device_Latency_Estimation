clear;
% load seattle_RTT;
% M=T(:,:,1);
% nmlz=max(M1(:));
% M=M1/nmlz;
M=rand(20);
[n1,n2]=size(M);
B=zeros(n1,n2);
beta=0.5;%采样率
m=round(n1*n2*beta);
si=randsample(n1*n2,m);
B(si)=1;
X=M.*B;
t=0.01;%误差容忍度
Fhat=IS_p(X,B,t);
Mhat=X+Fhat.*(1-B);
E=Mhat-M;
RSE=norm(E(:))/norm(M(:));
