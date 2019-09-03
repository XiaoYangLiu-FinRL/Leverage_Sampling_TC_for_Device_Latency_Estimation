function [ X4 ] = infmc( X,m)
%INFMC:MCE2E采样。 输入完整矩阵X和采样数m，返回采样并恢复后的矩阵X4.
%   限制：X为方阵
beta=0.4;
[N,~]=size(X);
m1=round(beta*m);
etr=round(0.5*N*log(N));
e=0.01;
sp=randsample(N*N,m1+etr);
B2=zeros(N*N,1);
B1=B2;
for nn=1:m1
    B1(sp(nn),1)=1;
end
B1=reshape(B1,N,N);
X1=B1.*X;
D1=sparse(X1);
[A1, ~, ~] = inexact_alm_mc(D1, 1e-3);
X2=A1.U*A1.V';
X2=X1+X2.*(1-B1);
for nn=1:m1+etr
    B2(sp(nn),1)=1;
end
B2=reshape(B2,N,N);
X2_=B2.*X;
D2=sparse(X2_);
[A2, ~, ~] = inexact_alm_mc(D2, 1e-3);
X3=A2.U*A2.V';
X3=X2_+X3.*(1-B2);
itr=2;
miu=0.01;
for ir=1:itr
   info=zeros(N,N);
   xita=info;
   alpha=0;
    for i=1:N
        for j=1:N
            info(i,j)=2*abs(X3(i,j)-X2(i,j))/abs(X3(i,j)+X2(i,j));
            if info(i,j)>miu
                xita(i,j)=1;
                alpha=alpha+1;
            end
        end
    end
%     info1=reshape(info,N1*N2,1);
    B2_temp=1-B2;
    info_remain=info.*B2_temp;
%      alpha=alpha/(N*N);
    [~,si]=sort(info_remain(:),'descend');
    new=round((m-m1-etr)/itr);
    for k=1:new
       B2(si(k))=1;
    end
    X2_=B2.*X;
    D2=sparse(X2_);
    [A2, ~, ~] = inexact_alm_mc(D2, 1e-3);
    X4=A2.U*A2.V';
    X4=X2_+X4.*(1-B2);
    e1=equal(X3,X4);
    if e1<=e
        break     
    end
    X2=X3;
    X3=X4;
end

end

