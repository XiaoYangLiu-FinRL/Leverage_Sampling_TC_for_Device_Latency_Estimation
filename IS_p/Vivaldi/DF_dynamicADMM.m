%% ��һ��Ƭ�ܸ˷����ϸߵĵ�����һ��Ƭ������������ֻռһ���ı���
clear;
load seattle_RTT;
l=T(:,:,200:250);
clear T;
T=l;
[n1,n2,n3]=size(T);
D=zeros(n1,n2,n3);
F=D;
normalize=max(T(:));
X=T/normalize;
r=3;
betanumber=0;
%  beta=0.6;%0.3Ч�����
rate_sample=0.1:0.1:0.9;
[~,rate_number]=size(rate_sample);
beta=0.3;
RSE=zeros(1,rate_number)  ;

for i=1:rate_number%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����
    m=round(0.1*i*n1*n2);
    B=zeros(n1,n2,n3);
    temp1=X(:,:,1);
    m1=ceil(m*beta);
    m2=m-m1;
    Omega1 = randsample(n1*n2,m1);
    M=zeros(n1,n2);
    M(Omega1) = temp1(Omega1);
    [probability] = leverage_scores(M,r);
    [~, si]=sort(probability(:),'descend'); % svΪ�������ֵ��siΪ��Ӧ��������
    [ remainer ] = fordifferent( si,Omega1 );
    Btemp=zeros(n1,n2);
    Btemp(Omega1)=1;
    Btemp(remainer(1:m2))=1;
    B(:,:,1)=Btemp;
    for k=2:n3
        M=X(:,:,k-1).*B(:,:,k-1);
        [probability] = leverage_scores(M,r);
        [~, si]=sort(probability(:),'descend');
        ri=randsample(n1*n2,n1*n2);
        Btemp=zeros(n1,n2);
        Btemp(si(1:m1))=1;
        remainer1=fordifferent(ri,si(1:m1));
        Btemp(remainer1(1:m2))=1;
        B(:,:,k)=Btemp;
        if k==n3
            break
        end
    end
    imcompleteRTT=X.*B;
    for d=1:n3
        temp=imcompleteRTT(:,:,d);
         D(:,:,d)= NCS_vivaldi_all(temp, 8, length(temp), 32, 0);
    end
    
    for ii=1:n1
        for jj=1:n2
            for kk=1:n3
                if D(ii,jj,kk)==0
                   F(ii,jj,kk)=0;
                else
                   F(ii,jj,kk)=imcompleteRTT(ii,jj,kk)/D(ii,jj,kk);
                end
            end
        end
    end
    
    
    % ADMM
    %M:��������    omega:����Ԫ�����     B:����Ԫ��Ϊ1,�ǳ���Ԫ��Ϊ0������
    alpha                  =        1.5;  %���ɳڲ���(alpha��ֵ����1.0��1.8֮��)
    maxItr                 =        50; % maximum iteration
    rho                    =        0.01; %�����������ղ���
    p                      =        0.5 ;
    addpath('tSVD','proxFunctions','solvers')
    myNorm  = 'tSVD_1' ;
    A=diag(sparse(double(B(:))));
    M1=A * F(:);        %A:n1n2n3*n1n2n3�ĳ���Ԫ��Ϊ1,�ǳ���Ԫ��Ϊ0       M1��n1n2n3*1�ĳ�������
    X1                     =        tensor_cpl_admm(A,M1,rho,alpha,[n1,n2,n3],maxItr,myNorm,0);
    X1                     =        reshape(X1,[n1,n2,n3])        ;
    X2=X1.*D;
    X2 = X2*normalize;

    X_dif                  =        X2-T                          ;
    RSE(i)               =        norm(X_dif(:))/(norm(T(:)));               %relative_square_error ���ƽ�����
end


%  open('alm_admm.fig');
figure;
k=0.1:0.1:0.9;
plot(k,RSE);
