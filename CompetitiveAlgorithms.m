%% 对比算法，包括Random sampling+admm, Random sampling+alm, MCE2E+alm.
clc;
clear;
tic
%  T=fpower_law(2);%产生powerlaw分布的低秩数据
% [n1,n2,n3]=size(T1);
%加噪
% signal_power = 1/length(T1(:))*sum(T1(:).*T1(:));
% T2=1/sqrt(signal_power) * T1;

% T2=T1/norm(T1(:));
% R=normrnd(0,0.002,n1,n2,n3);
% T=T2+R;

% alpha=0.6;
% idx=randsample(n1*n2*n3,round(alpha*n1*n2*n3));
% long=size(idx);
% for f=1:long
%    T(idx(f)) =T(idx(f))+T(idx(f))*0.2*rand;
% end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%
% load PlanetLab_RTT;
% o=T(1:120,1:120,1:3);
% clear T;
% T=o;
% clear;
% Nway = [50,50,10]; % dimension of tensor
% r = 5;
% coreNway = [r,r,r];
% % 
% % % randomly generate core tensor
% rr=rand(coreNway);
% G = tensor(rr);
% A = cell(1,ndims(G));
% % randomly generate factor matrices
% for i = 1:ndims(G)                   % ndims(G) 张量阶数
% a=randn(Nway(i),r);
% b= orth(a);
% A{i} =b*diag((1:r).^(-0.5));
% end
% % generate tensor
% X2= full(ttensor(G,A));
% T=X2.data;
% clearvars -EXCEPT T;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% load PowerLaw_50_50_10e2.mat;
load Seattle_RTT.mat
q=T(:,:,1:10);%取十帧为例进行实验
clearvars T;
T=q;
[n1,n2,n3]=size(T);
Y=T(:);
% signal_power = 1/length(Y)*sum(Y.*Y);
% signal_power=norm(Y,'fro')^2;
% normalize              =  sqrt(signal_power);
normalize=max((Y));
X=T/normalize;
% X=(X+(normrnd(0,0.1,[n1,n2,n3])));%加噪
r=5;
rho  =  0.01;
RSE1=zeros(1,9);
RSE2=zeros(1,9);
RSE3=zeros(1,9);
RSE4=zeros(1,9);
RSE1=zeros(1,9);
RSE2sum=RSE1;
RSE3sum=RSE1;
RSE4sum=RSE1;
sample_rate=0.1:0.1:0.9;
circletimes=10;%循环多次取均值
for ct=1:circletimes
    for k=1:9%采样率0.1-0.9
        
%         mm1=3;         %随机采样的片数
%         B2=zeros(n1,n2,n3);
%         num=round(n1*n2*k*0.1);         %每一片的采样数
%         for m=1:mm1
%             sp=randsample(n1*n2,num);
%             B1=zeros(n1,n2);
%             B1(sp)=1;
%             B2(:,:,m)=B1;
%         end
%         for p=mm1+1:n3
%             X1=X.*B2;
%             [UU,SS,VV]=reduce_tsvd(X1,r);
%             U_h_slide={1,n1};
%             V_l_slide={1,n2};
%             u=zeros(1,n1);   %存储行相关系数,左奇异值的正切片
%             v=zeros(1,n2);   %存储列相关系数，右奇异值的侧切片
%             %水平切片 horizon slide
%             for i=1:n1
%                 U_h_slide{1,i}=reshape(UU(i,:,:),r,n3);
%                 u(1,i)=norm(U_h_slide{1,i}(:))^2;
%                 %      u(1,i)=norm(U_h_slide{1,i},2)^2;
%             end
%             %侧切片   lateral slide
%             for j=1:n2
%                 V_l_slide{1,j}=reshape(VV(j,:,:),r,n3);
%                 v(1,j)=norm(V_l_slide{1,j}(:))^2;
%                 %     v(1,j)=norm(V_l_slide{1,j},2)^2;
%             end
%             p1=zeros(n1,n2);    %杠杆分数
%             
%             for i=1:n1
%                 for j=1:n2
%                     p1(i,j)=(u(1,i)*r/n1+v(1,j)*r/n2-(u(1,i)*r/n1)*(v(1,j)*r/n2))*(log(n1+n2)^2);
%                 end
%             end
%             p1_array=reshape(p1,[1,n1*n2]);
%             [p1_descend,Idx]=sort(p1_array,'descend');
%             %         remainer=fordifferent(Idx,sp);
%             
%             B3=zeros(1,n1*n2);
%             B3(Idx(1:num))=1;
%             %         for ii=1:round(mm*(1-beta))
%             %             B3(1,remainer(ii))=1;
%             %         end
%             B4=reshape(B3,[n1,n2]);
%             %             B5=zeros(n1,n2,n3);
%             %             for kk=mm+1:n3
%             %                 B5(:,:,kk)=B4; %杠杆采样索引
%             %             end
%             %             BB=B2+B5;%总索引
%             B2(:,:,p)=B4;
%         end
%         BB=B2;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %% ADMM
%         % M:抽样张量    omega:抽样元素序号     B:抽样元素为1,非抽样元素为0的张量
%         alpha                  =        1.5;  %过松弛参数(alpha的值是在1.0和1.8之间)
%         maxItr                 =        50; % maximum iteration
%         rho                    =        0.01; %增广拉格朗日参数
%         p                      =        0.5 ;
%         addpath('tSVD','proxFunctions','solvers')
%         myNorm  = 'tSVD_1' ;
%         A=diag(sparse(double(BB(:))));
%         M1=A * X(:);        %A:n1n2n3*n1n2n3的抽样元素为1,非抽样元素为0       M1：n1n2n3*1的抽样向量
%         X1                     =        tensor_cpl_admm(A,M1,rho,alpha,[n1,n2,n3],maxItr,myNorm,0);
%         X1 = X1*normalize;
%         Xleverage                     =        reshape(X1,[n1,n2,n3])        ;
%         X_dif                  =        Xleverage-T                          ;
%         RSE1(1,k)               =        norm(X_dif(:))/(norm(T(:)));
        mm=round(0.1*n1*n2*k);         %每片样本大小设置
        %% 随机采样+admm
        rB2=zeros(n1,n2,n3);
        
        for kk=1:n3
            rB2_temp=zeros(n1,n2);
            sp=randsample(n1*n2,mm);
            rB2_temp(sp)=1;
            rB2(:,:,kk)=rB2_temp;
        end
        % ADMM
        %M:抽样张量    omega:抽样元素序号     B:抽样元素为1,非抽样元素为0的张量
        alpha                  =        1.5;  %过松弛参数(alpha的值是在1.0和1.8之间)
        maxItr                 =        50; % maximum iteration
        rho                    =        0.01; %增广拉格朗日参数
        p                      =        0.5 ;
        addpath('tSVD','proxFunctions','solvers')
        myNorm  = 'tSVD_1' ;
        A=diag(sparse(double(rB2(:))));%X是原始完整数据，rB2是索引。
        Xtemp1=X.*rB2; %Xtemp1是采样张量
        M1=A * Xtemp1(:);        %A:n1n2n3*n1n2n3的抽样元素为1,非抽样元素为0       M1：n1n2n3*1的抽样向量
        X1                     =        tensor_cpl_admm(A,M1,rho,alpha,[n1,n2,n3],maxItr,myNorm,0);
        X1 = X1*normalize;
        Xrandtube                     =        reshape(X1,[n1,n2,n3])        ;
        X_dif                  =        Xrandtube-T                          ;
        RSE2(1,k)               =        norm(X_dif(:))/(norm(T(:)));
        %%  MCE2E+alm
        T_INF=zeros(n1,n2,n3);
        for i=1:n3
            Xtemp=X(:,:,i);
            T_INF(:,:,i)=infmc(Xtemp,mm);
        end
        X_dif                  =        T_INF*normalize-T                          ;
        RSE3(1,k)               =        norm(X_dif(:))/(norm(T(:)));
        %% 随机采样+ALM恢复
        Xrand=zeros(n1,n2,n3);
        for ii=1:n3
            Xtemp=X(:,:,ii);
            sp_rand=randsample(n1*n2,mm);
            B_rand=zeros(n1,n2);
            B_rand(sp_rand)=1;
            Xr=B_rand.*Xtemp;
            D1=sparse(Xr);
            [A1, ~, ~] = inexact_alm_mc(D1, 1e-3);
            X2=A1.U*A1.V';
            X2=Xr+X2.*(1-B_rand);
            Xrand(:,:,ii)=X2;
        end
        X_dif                  =        Xrand*normalize-T                          ;
        RSE4(1,k)               =        norm(X_dif(:))/(norm(T(:)));
        
        
    end
    
%     RSE1sum=RSE1sum+RSE1;
    RSE2sum=RSE2sum+RSE2;
    RSE3sum=RSE3sum+RSE3;
    RSE4sum=RSE4sum+RSE4;
    
end

%% 画图

figure;
% plot(sample_rate,RSE1sum/circletimes,'r');
% hold on;
plot(sample_rate,RSE2sum/circletimes,'b');
hold on;
plot(sample_rate,RSE3sum/circletimes,'k');
hold on;
plot(sample_rate,RSE4sum/circletimes,'r--');
% hold on;
% load RSE_PlanetLab_RTT_120_120_3;
% plot(sample_rate,RSE,'k--*');
xlabel('sample rate');
ylabel('RSE');
legend('Random sampling+ADMM','MCE2E+ALM','Random sampling+ALM');
% legend('Location','southeast')
grid on;
toc





