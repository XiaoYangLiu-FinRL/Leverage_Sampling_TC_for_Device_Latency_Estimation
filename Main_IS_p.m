%% IS-p算法，需要安装Mosek来解凸优化。默认取p=2
clear;
% cvx_setup;
tic
load PlanetLab_RTT.mat;
M=T(:,:,1:3);
[n1,n2,n3]=size(M);
default_dimension=8;

RSE=zeros(1,9);
for beta1=3:9
    beta=0.1*beta1;%采样率
    B=zeros(n1,n2,n3);
    D=zeros(n1,n2,n3);
    m=round(n1*n2*n3*beta);
    si=randsample(n1*n2*n3,m);
    B(si)=1;
    X=M.*B;
    for k=1:n3
        Xtemp=X(:,:,k);
        D1 = NCS_vivaldi_all(Xtemp, default_dimension, length(Xtemp), 32, 0);
        D(:,:,k)=D1;
        
    end
    %     Dtemp=D;
    % %     idx=find(D==0);
    %     Dtemp(D==0)=0.001;
    F=X./D;
    F1 = fillmissing(F,'constant',0);
    F2=F1;
    F2(F1==inf)=0;
    t=0.1;%误差容忍度
    Fhat=TC_IS_p(F2,B,t);
    Mhat=Fhat.*D;
    % Mhat=abs(Mhat);
    M2=X+Mhat.*(1-B);
    E=M2-M;
    RSE(1,beta1)=norm(E(:))/norm(M(:));

end
toc