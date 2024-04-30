close all;clear;
%多个展开点检测性能对比曲线
tic
pf=1e-4;
K=[32 29 26 23];%reference number
%K=[32 32 32 32];
L=[16 16 16 16];%snapshot
site_num=length(K);
M=K+1-L;
rou=(K+2-L)./(K+1);%损失因子均值 
w_modify=(K+1)./(K+1-L);
trial_num0=1e6;
trial_num1=1e6; 
SNR=0:1:20;
lamda_snr=10.^(SNR/10);
snr_num=length(lamda_snr);

%load('cof_opt1111.mat') %10621
%load('cof_opt1113.mat') %1711
%load('cof_opt1112.mat') %11 1415
load('cof_opt1114.mat') %1215

cof_opt1111=cof_opt1114; %更改一下参数名
 
% lamda_snrtimes=[10 6 2 1]; %局部信噪比的比例
 lamda_snrtimes=[1 2 1 5]; %局部信噪比的比例
%  lamda_snrtimes=[1 7 1 1]; %局部信噪比的比例
% lamda_snrtimes=[1 1 14 15]; %局部信噪比的比例 %用公式求门限在某些展开点是不正确的


snr_matrix=[];
for i=1:length(lamda_snrtimes)
    lamda_snr_local=lamda_snrtimes(i)*lamda_snr/sum(lamda_snrtimes);
    beta_snr_local=lamda_snr_local/rou(i);
    snr_matrix=[snr_matrix beta_snr_local'];
end


cof_opt=zeros(1,site_num);
for i=1:site_num
    cof_opt(i)=cof_opt1111(i);
end

load('MultiExpandValue.mat')
eta04=eta04(1,:);
expand_num=size(eta02,1);

pd_expand_r1=zeros(1,snr_num);pd_expand_r2=zeros(expand_num,snr_num);
pd_expand_r3=zeros(1,snr_num);pd_expand_opt=zeros(1,snr_num);
pd_standardZ=zeros(1,snr_num);pd_modifyZ=zeros(1,snr_num);

gate_r2=zeros(1,expand_num);
gatee_r2=zeros(1,expand_num);


cof_r1=zeros(1,site_num);cof_r2=zeros(expand_num,site_num);
cof_r3=zeros(1,site_num);

%%
for h=1:expand_num
    t0_r1=0;t0_r2=0;t0_r3=0;%泰勒展开
    t0_opt=0;
    t0_standardZ=0;t0_modifyZ=0;
    temp1=zeros(1,site_num);
    temp2=zeros(1,site_num);
    temp3=zeros(1,site_num);
    cof_eta=zeros(1,site_num);
    cof_z=zeros(1,site_num);
    for i=1:site_num
        %求各检测器权值
        temp1(i)=hypergeom(K(i)+3-L(i),2,snr_matrixE(i).*rou(i).*eta01(i))....
            /hypergeom(K(i)+2-L(i),1,snr_matrixE(i).*rou(i).*eta01(i)); %均值
        temp2(i)=hypergeom(K(i)+3-L(i),2,snr_matrixE(i).*rou(i).*eta02(h,i))....
            /hypergeom(K(i)+2-L(i),1,snr_matrixE(i).*rou(i).*eta02(h,i)); %2H1分位点/3H0分位点
        temp3(i)=hypergeom(K(i)+3-L(i),2,snr_matrixE(i).*rou(i).*eta04(i))....
            /hypergeom(K(i)+2-L(i),1,snr_matrixE(i).*rou(i).*eta04(i)); %中位数
        %         cof_eta(i)=(K(i)+2-L(i)).*snr_matrix(k,i).*rou(i).*temp1(i);%eta
        %         cof_z(i)=cof_eta(i).*(1-eta01(i)).^2;
        cof_r1_temp(i)=rou(i).^2.*snr_matrixE(i).*(1-eta01(i)).*temp1(i);%r1
        cof_r2_temp(h,i)=rou(i).^2.*snr_matrixE(i).*(1-eta02(h,i)).*temp2(i);%r2
        cof_r3(i)=rou(i).^2.*snr_matrixE(i).*(1-eta04(i)).*temp3(i);%r1
%         %cof_x1(i)=cof_z(i)/(K(i)+1)*(1-eta01(i))^(K(i));
        %H0下统计量MC
        randsig=(randn(M(i)+1,trial_num0)+1i*randn(M(i)+1,trial_num0))/sqrt(2); %CN(0,I),噪声方差为1
        r=(K(i)+1)*log(1+abs(randsig(1,:)).^2./sum(abs(randsig(2:end,:)).^2,1));%r
        
%         cof_r1=round(100*(cof_r1/cof_r1(1)))/100;
%         cof_r2=round(100*(cof_r2./cof_r2(:,1)))/100;
        cof_r1(i)=(cof_r1_temp(i)/cof_r1_temp(1));
        cof_r2(h,i)=cof_r2_temp(h,i)/cof_r2_temp(h,1);
        
%         t0_r1=t0_r1+cof_r1(i)*r;%r-expand
%         t0_r2=t0_r2+cof_r2(h,i)*r;%r-expand
%         t0_r3=t0_r3+cof_r3(i)*r;
%         t0_opt=t0_opt+cof_opt(i)*r;
%         
%         t0_standardZ=t0_standardZ+r;%standard
%         t0_modifyZ=t0_modifyZ+2/w_modify(i)*r;%modify
    end
%     t0_r1=sort(t0_r1);
%     t0_r2=sort(t0_r2);
%     t0_opt=sort(t0_opt);
%     t0_r3=sort(t0_r3);
% 
%     t0_standardZ=sort(t0_standardZ);t0_modifyZ=sort(t0_modifyZ);
%     %%
%     gate_r1=t0_r1(trial_num0-trial_num0*pf);
%     gate_r2(h)=t0_r2(trial_num0-trial_num0*pf);
%     gate_r3=t0_r3(trial_num0-trial_num0*pf);
%     gate_opt=t0_opt(trial_num0-trial_num0*pf);
% 
%     gate_standardZ=t0_standardZ(trial_num0-trial_num0*pf);
%     gate_modifyZ=t0_modifyZ(trial_num0-trial_num0*pf);
            % 用公式求门限
            gatee_r1=wchigate(2*ones(site_num,1),(w_modify.*cof_r1)/2,pf,28);
            gatee_r2(h)=wchigate(2*ones(site_num,1),(w_modify.*cof_r2(h,:))/2,pf,28);
            gatee_r3=wchigate(2*ones(site_num,1),(w_modify.*cof_r3)/2,pf,28);
            gatee_opt=wchigate(2*ones(site_num,1),(w_modify.*cof_opt)/2,pf,28);

            
            gatee_standardZ=wchigate(2*ones(site_num,1),(w_modify)/2,pf,28);
            gatee_modifyZ=chi2inv(1-pf,2*site_num);
end
cof_r2aaa1=cof_r2./cof_r2(:,1);
cof_r1aaa1=cof_r1/cof_r1(1);
ddr2=cof_r2-cof_r1;
ddcofr2=sum(ddr2.^2,2);
% for k=1:snr_num %%这两行只能看看门限计算是否正确 画检测曲线有问题（系数c）
%k/snr_num*100
%计算检测概率
for h=1:expand_num
    h/expand_num
    for k=1:snr_num
        t1_r1=0;t1_r2=0;t1_r3=0;t1_opt=0;
        t1_standardZ=0;t1_modifyZ=0;
        for i=1:site_num
            randsig=(randn(M(i)+1,1,trial_num1)+1i*randn(M(i)+1,1,trial_num1))/sqrt(2);
            loss_factor=betarnd(K(i)-L(i)+2,L(i)-1,1,trial_num1);
            r=(K(i)+1)*log(1+abs(randsig(1,:)+sqrt(loss_factor)*sqrt(snr_matrix(k,i))).^2./sum(abs(randsig(2:end,:)).^2));
            t1_r1=t1_r1+cof_r1(i)*r;
            t1_r2=t1_r2+cof_r2(h,i)*r;
            t1_r3=t1_r3+cof_r3(i)*r;
            t1_opt=t1_opt+cof_opt(i)*r;
            
            t1_standardZ=t1_standardZ+r;
            t1_modifyZ=t1_modifyZ+2/w_modify(i)*r;
        end
%         pd_expand_r1(k)=sum(t1_r1>gate_r1);
%         pd_expand_r2(h,k)=sum(t1_r2>gate_r2(h));
%         pd_standardZ(k)=sum(t1_standardZ>gate_standardZ);
%         pd_modifyZ(k)=sum(t1_modifyZ>gate_modifyZ);

        pd_expand_r1(k)=sum(t1_r1>gatee_r1);
        pd_expand_r2(h,k)=sum(t1_r2>gatee_r2(h));
        pd_expand_r3(k)=sum(t1_r3>gatee_r3);
        pd_expand_opt(k)=sum(t1_opt>gatee_opt);
        
        pd_standardZ(k)=sum(t1_standardZ>gatee_standardZ);
        pd_modifyZ(k)=sum(t1_modifyZ>gatee_modifyZ);
    end
end
pd_expand_opt=pd_expand_opt/trial_num1*100;
pd_expand_r1=pd_expand_r1/trial_num1*100;
pd_expand_r2=pd_expand_r2/trial_num1*100;
pd_expand_r3=pd_expand_r3/trial_num1*100;
pd_standardZ=pd_standardZ/trial_num1*100;
pd_modifyZ=pd_modifyZ/trial_num1*100;
figure; 
linesty={'-*','-+','-o','-d','-s','-<','-x','-v'};
plot(SNR,pd_expand_opt,'-h','linewidth',1.5,'markersize',8);hold on;
for h=1:expand_num
    plot(SNR,pd_expand_r2(h,:),linesty{h},'linewidth',1.5,'markersize',8);hold on;
end
plot(SNR,pd_expand_r1,'->','linewidth',1.5,'markersize',8);hold on;
%plot(SNR,pd_expand_r3,'-v','linewidth',1.5,'markersize',8);hold on;
% plot(SNR,pd_standardZ,'-p','linewidth',1.5,'markersize',8);hold on;
% plot(SNR,pd_modifyZ,'-x','linewidth',1.5,'markersize',8);
grid on;
xlabel('TSNR (dB)');ylabel('Probability of detection (%)')

% % legend('\alpha_{0}=0.01%','\alpha_{0}=0.1%','\alpha_{0}=1%','\alpha_{0}=10%'...
% %      ,'\alpha_{0}=20%','\alpha_{0}=50%','\alpha_{0}=90%','Mean:\lambda_0=9dB');
legend('Opt','\alpha_{1}=1%','\alpha_{1}=10%','\alpha_{1}=30%','\alpha_{1}=40%'...
     ,'\alpha_{1}=50%','\alpha_{1}=70%','\alpha_{1}=90%','Mean:\lambda_0=9dB');
set(gcf,'color',[1,1,1]);
set(gca,'Fontname','Times New Roman','FontSize',13);

toc