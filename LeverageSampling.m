%% 杠杆采样。
clc;
clear;
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load seattle_RTT.mat;
o=T(:,:,1:10);
clear T;
T=o;
%%%%%%%%%%%%%%%%%%%%%%
% T=fpower_law(2);
% load PowerLaw_50_50_10.mat;
% o=abs(T);
% clear T;
% T=o;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load PlanetLab_RTT;
% o=T(1:120,1:120,1:10);
% clear T;
% T=o;
%%
r=2;
beta=0.1;               %随机采样的片数所占比例
rho  =  0.01;
[n1,n2,n3]=size(T);
normalize              =        max(abs(T(:)))                     ;
RSE1=zeros(1,9);
RSE1sum=RSE1;
sample_rate=0.1:0.1:0.9;
circletimes=5;
for ct=1:circletimes
    for k=1:9
        X=T/normalize;
        mm=round(n3*beta);         %随机采样的片数
        B2=zeros(n1,n2,n3);
        num=round(n1*n2*k*0.1);         %每一片的采样数
        for m=1:mm
            sp=randsample(n1*n2,num);
            B1=zeros(n1,n2);
            B1(sp)=1;
            B2(:,:,m)=B1;
        end
        for p=mm+1:n3
            X1=X.*B2;
            [UU,SS,VV]=reduce_tsvd(X1,r);
            U_h_slide={1,n1};
            V_l_slide={1,n2};
            u=zeros(1,n1);   %存储行相关系数,左奇异值的正切片
            v=zeros(1,n2);   %存储列相关系数，右奇异值的侧切片
            %水平切片 horizon slide
            for i=1:n1
                U_h_slide{1,i}=reshape(UU(i,:,:),r,n3);
                u(1,i)=norm(U_h_slide{1,i}(:))^2;
                %      u(1,i)=norm(U_h_slide{1,i},2)^2;
            end
            %侧切片   lateral slide
            for j=1:n2
                V_l_slide{1,j}=reshape(VV(j,:,:),r,n3);
                v(1,j)=norm(V_l_slide{1,j}(:))^2;
                %     v(1,j)=norm(V_l_slide{1,j},2)^2;
            end
            p1=zeros(n1,n2);    %杠杆分数
            
            for i=1:n1
                for j=1:n2
                    p1(i,j)=(u(1,i)*r/n1+v(1,j)*r/n2-(u(1,i)*r/n1)*(v(1,j)*r/n2))*(log(n1+n2)^2);
                end
            end
            p1_array=reshape(p1,[1,n1*n2]);
            [p1_descend,Idx]=sort(p1_array,'descend');
            %         remainer=fordifferent(Idx,sp);
            
            B3=zeros(1,n1*n2);
            B3(Idx(1:num))=1;
            %         for ii=1:round(mm*(1-beta))
            %             B3(1,remainer(ii))=1;
            %         end
            B4=reshape(B3,[n1,n2]);
            %             B5=zeros(n1,n2,n3);
            %             for kk=mm+1:n3
            %                 B5(:,:,kk)=B4; %杠杆采样索引
            %             end
            %             BB=B2+B5;%总索引
            B2(:,:,p)=B4;
        end
        BB=B2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ADMM
        % M:抽样张量    omega:抽样元素序号     B:抽样元素为1,非抽样元素为0的张量
        alpha                  =        1.5;  %过松弛参数(alpha的值是在1.0和1.8之间)
        maxItr                 =        50; % maximum iteration
        rho                    =        0.01; %增广拉格朗日参数
        p                      =        0.5 ;
        addpath('tSVD','proxFunctions','solvers')
        myNorm  = 'tSVD_1' ;
        A=diag(sparse(double(BB(:))));
        M1=A * X(:);        %A:n1n2n3*n1n2n3的抽样元素为1,非抽样元素为0       M1：n1n2n3*1的抽样向量
        X1                     =        tensor_cpl_admm(A,M1,rho,alpha,[n1,n2,n3],maxItr,myNorm,0);
        X1 = X1*normalize;
        Xleverage                     =        reshape(X1,[n1,n2,n3])        ;
        X_dif                  =        Xleverage-T                          ;
        RSE1(1,k)               =        norm(X_dif(:))/(norm(T(:)));
        
    end
    
    RSE1sum=RSE1sum+RSE1;
end
figure;
plot(sample_rate,RSE1sum/circletimes,'r');

xlabel('sample rate');
ylabel('RSE');
legend('Leverage sampling+ADMM');
% legend('Location','southeast')
grid on;
toc