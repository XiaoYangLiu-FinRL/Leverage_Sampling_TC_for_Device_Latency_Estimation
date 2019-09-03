%张量每片都做矩阵score采样
clc;
clear;
% T=importdata('SeattleData_1');
load seattle_RTT.mat;
p=T(:,:,1);
clear T;
T=p;
%基本参数设置
r=3; %r=10还行
% n1=100;
% n2=100;
[n1,n2]=size(T);
% error_bound=0.01; %一次试验中，当误差小于多少，认为恢复正确
Beta=0.5; %0.5还行
% alph=0.5;
tt=10;  %运行多次求均值
% L=8;
% %采样率改变
rate_sample=0.1:0.1:0.9;                  %采样率
[~,rate_number]=size(rate_sample);
default_dimension=8;
T1_rse=zeros(tt,rate_number);
for i=1:tt
    for l=1:rate_number
        T_fnorm=norm(T(:));
        m=round(0.1*l*n1*n2);       %总的样本数
        %观测集合
        [M1,B1,Omega1]=M_Npass_sampling(m,T,r,4,Beta);
        %     [M2,B2,Omega2]=M_uniform(m,X);
        F=zeros(size(M1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         D = NCS_vivaldi_all(M1, default_dimension, length(M1), 32, 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 [out_all, in_all, fpre_newhost, fpre_flashcrowd] = NCS_phoenix(M1, default_dimension, length(M1), 32, 5, 0, []);        
%         D = out_all * in_all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                [out_all, in_all, fperr] = NCS_DMF(M1, default_dimension, 32, 100, 50, 0);
%         D = out_all * in_all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = length(T);
%         tmp = randperm(N);         
%         landmarks = tmp(1:32); % choose random L nodes as landmarks
%         hosts = tmp(32+1:N);  % the rest of them are orinary hosts
%         D_landmark = T(landmarks, landmarks);
%         D_host2landmark = T(hosts, landmarks);
%                [out_l, in_l, out_h, in_h] = NCS_IDES_all(D_landmark, D_host2landmark, default_dimension, 0);
%         D = out_h*in_h;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load vivaldi.mat;
D=predicted_matrix;
        for j=1:n1
            for k=1:n2
                if D(j,k)==0
                    F(j,k)=0;
                else
                    F(j,k)=M1(j,k)/D(j,k);
                end
            end
        end
        normalize              =        max(abs(F(:)))    ;
        X=F/normalize;
        D1=sparse(X);
        [A1, iter1, svp1] = inexact_alm_mc(D1, 1e-4);
        X1=A1.U*A1.V';
        X1=X+X1.*(1-B1);
        Fhat=X1*normalize;
        Mhat=Fhat.*D;
        E=Mhat-T;
        T1_rse(i,l)=norm(E(:))/norm(T(:));
    end
    
end
r1=zeros(1,rate_number);

for k=1:tt
    r1=r1+T1_rse(k,:);
end
rse1=r1/tt;

%画图
%采样率和正确恢复比例的关系
figure(1);
plot(rate_sample,rse1,'r');

xlabel('rate sample');
ylabel('RSE');
legend('DF-ALM');






