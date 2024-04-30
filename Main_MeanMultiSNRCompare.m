close all;clear;
%多个展开点检测性能对比曲线
%在不同的总信噪比下
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
trial_num1=2e6;
SNR=0:1:20;
lamda_snr=10.^(SNR/10);
snr_num=length(lamda_snr);
% 
lamda_snrtimes=[1 2 1 5;
               10 6 2 1;
               1 7 1 1;
               1 1 14 15]; %局部信噪比的比例

% lamda_snrtimes=[1 2 1 5]; %局部信噪比的比例
% lamda_snrtimes(m,:)=[1 7 1 1]; %局部信噪比的比例
% lamda_snrtimes(m,:)=[1 1 14 15]; %局部信噪比的比例
fig_num=size(lamda_snrtimes,1);
for m=1:fig_num
snr_matrix=[];
for i=1:length(lamda_snrtimes(m,:))
    lamda_snr_local=lamda_snrtimes(m,i)*lamda_snr/sum(lamda_snrtimes(m,:));
    beta_snr_local=lamda_snr_local/rou(i);
    snr_matrix=[snr_matrix beta_snr_local'];
end

if m==1
    load('MultiMeanValue1215.mat')
elseif m==2
    load('MultiMeanValue10621.mat')
elseif m==3
    load('MultiMeanValue1711.mat')
else
    load('MultiMeanValue111415.mat')
end


mean_expand_all=eta01;
choice_index=[6 8 10 12 14]; %对应 5  7 9 11 13 dB
%choice_index=10;
%mean_expand_choice=eta01(:,[1 5 9 13 17 21]);
mean_expand_choice=eta01(:,choice_index);
%mean_expand_choice=eta01;
mean_expand_choice_num=size(mean_expand_choice,2);

pd_expand_r1=zeros(mean_expand_choice_num,snr_num);
pd_standardZ=zeros(1,snr_num);pd_modifyZ=zeros(1,snr_num);

gate_r2=zeros(1,mean_expand_choice_num);
gatee_r2=zeros(1,mean_expand_choice_num);


cof_r1=zeros(site_num,mean_expand_choice_num);
gatee_r1=zeros(1,mean_expand_choice_num);
gate_r1=zeros(1,mean_expand_choice_num);
%%
for h=1:mean_expand_choice_num
    t0_r1=0;t0_r2=0;t0_r3=0;%泰勒展开
    t0_standardZ=0;t0_modifyZ=0;
    temp1=zeros(1,site_num);
    for i=1:site_num
        %求各检测器权值
        temp1(i)=hypergeom(K(i)+3-L(i),2,snr_matrix(choice_index(h),i).*rou(i).*mean_expand_choice(i,h))....
            /hypergeom(K(i)+2-L(i),1,snr_matrix(choice_index(h),i).*rou(i).*mean_expand_choice(i,h)); %均值
        cof_r1(i,h)=rou(i).^2.*snr_matrix(choice_index(h),i).*(1-mean_expand_choice(i,h)).*temp1(i);%r1
        %H0下统计量MC
        randsig=(randn(M(i)+1,trial_num0)+1i*randn(M(i)+1,trial_num0))/sqrt(2); %CN(0,I),噪声方差为1
        r=(K(i)+1)*log(1+abs(randsig(1,:)).^2./sum(abs(randsig(2:end,:)).^2,1));%r        
        t0_r1=t0_r1+cof_r1(i,h)*r;%r-expand         
        t0_standardZ=t0_standardZ+r;%standard
        t0_modifyZ=t0_modifyZ+2/w_modify(i)*r;%modify
    end
    t0_r1=sort(t0_r1);
    t0_standardZ=sort(t0_standardZ);t0_modifyZ=sort(t0_modifyZ);
    %%
    gate_r1(h)=t0_r1(trial_num0-trial_num0*pf);
    gate_standardZ=t0_standardZ(trial_num0-trial_num0*pf);
    gate_modifyZ=t0_modifyZ(trial_num0-trial_num0*pf);
            % 用公式求门限
            gatee_r1(h)=wchigate(2*ones(site_num,1),(w_modify.*(cof_r1(:,h))')/2,pf,28);
            gatee_standardZ=wchigate(2*ones(site_num,1),(w_modify)/2,pf,28);
            gatee_modifyZ=chi2inv(1-pf,2*site_num);
end
% for k=1:snr_num %%这两行只能看看门限计算是否正确 画检测曲线有问题（系数c）
%k/snr_num*100
%计算检测概率
for h=1:mean_expand_choice_num
    h/mean_expand_choice_num
    for k=1:snr_num
        t1_r1=0; 
        t1_standardZ=0;t1_modifyZ=0;
        for i=1:site_num
            randsig=(randn(M(i)+1,1,trial_num1)+1i*randn(M(i)+1,1,trial_num1))/sqrt(2);
            loss_factor=betarnd(K(i)-L(i)+2,L(i)-1,1,trial_num1);
            r=(K(i)+1)*log(1+abs(randsig(1,:)+sqrt(loss_factor)*sqrt(snr_matrix(k,i))).^2./sum(abs(randsig(2:end,:)).^2));
            t1_r1=t1_r1+cof_r1(i,h)*r;
            t1_standardZ=t1_standardZ+r;
            t1_modifyZ=t1_modifyZ+2/w_modify(i)*r;
        end
%         pd_expand_r1(k)=sum(t1_r1>gate_r1);
%         pd_expand_r2(h,k)=sum(t1_r2>gate_r2(h));
%         pd_standardZ(k)=sum(t1_standardZ>gate_standardZ);
%         pd_modifyZ(k)=sum(t1_modifyZ>gate_modifyZ);

        pd_expand_r1(h,k)=sum(t1_r1>gatee_r1(h));
        pd_standardZ(k)=sum(t1_standardZ>gatee_standardZ);
        pd_modifyZ(k)=sum(t1_modifyZ>gatee_modifyZ);
    end
end
pd_expand_r1=pd_expand_r1/trial_num1*100;
pd_standardZ=pd_standardZ/trial_num1*100;
pd_modifyZ=pd_modifyZ/trial_num1*100;

figure;
linesty={'-*','-+','-o','-d','-s','-<','-x',':','-.'};
for h=1:mean_expand_choice_num
    plot(SNR,pd_expand_r1(h,:),linesty{h},'linewidth',1.5,'markersize',8);hold on;
end
%plot(SNR,pd_expand_r3,'-v','linewidth',1.5,'markersize',8);hold on;
% plot(SNR,pd_standardZ,'-p','linewidth',1.5,'markersize',8);hold on;
% plot(SNR,pd_modifyZ,'-x','linewidth',1.5,'markersize',8);
grid on;
xlabel('TSNR (dB)');ylabel('Probability of detection (%)')
%legend('5dB','9dB','13dB','17dB','21dB','25dB');
legend('$\hat{\lambda}_0$=5dB','$\hat{\lambda}_0$=7dB','$\hat{\lambda}_0$=9dB','$\hat{\lambda}_0$=11dB','$\hat{\lambda}_0$=13dB','Interpreter','latex');
%legend('SW','Standard GLRT','MGLRT');
set(gcf,'color',[1,1,1]);
set(gca,'Fontname','Times New Roman','FontSize',13);
end
toc