clear;
load seattle_RTT;
M=T(:,:,1);
% nmlz=max(M1(:));
% M=M1/nmlz;
% M=rand(20);
[n1,n2]=size(M);
B=zeros(n1,n2);
C=B;

beta=0.5;
s=round(n2*beta);
C1=zeros(n1,s);
X=M;
for i=1:s
    for j=1:n2
        temp=X(:,j);
        c(j)=(norm(temp(:),2))^2;
    end
    %     pr=c/sum(c);
    [~,si]=sort(c);
    if i~=1
        ti=si;
        ti=fordifferent(si,used) ;
        clear si;
        si=ti;
    end
    
    C(:,si(1))=M(:,si(1));
    B(:,si(1))=1;
    C1(:,i)=M(:,si(1));
    used(i)=si(1);
    X=X-C1*pinv(C1)*X;
    
    
    
    
end
D2=sparse(C);
[A2, iter2, svp2] = inexact_alm_mc(D2, 1e-4);
X2=A2.U*A2.V';
X2=C+X2.*(1-B);
E=X2-M;
rse=norm(E(:))/(norm(M(:)));
















