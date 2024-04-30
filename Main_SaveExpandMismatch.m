close all;clear;
%保存多个展开点
%均值以及H1下分位点
%H0下分位点
tic
pf=1e-4;
K=[32 29 26 23];%reference number
%K=[32 32 32 32];%reference number
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
AssumeAllSnr=10;
sum_snr=lamda_snr(AssumeAllSnr);

lamda_snrtimes=[10 6 2 1;1 1 1 1;1 2 3 4]; %局部信噪比的比例 %
%lamda_snrtimes=[1 1 14 15;1 1 1 1;4 3 2 1];
ratio_mismatch_num=size(lamda_snrtimes,1);
mismatch_num=ratio_mismatch_num+2; %总信噪比估计有差异 %额外增加两个总信噪比有偏差的曲线
lamda_snr_local=zeros(mismatch_num,site_num); %估计信噪比
%总信噪比已知，局部信噪比比值未知
for i=1:ratio_mismatch_num
    lamda_snr_local(i,:)=lamda_snrtimes(i,:)*sum_snr/sum(lamda_snrtimes(i,:));
    beta_snr_local(i,:)=lamda_snr_local(i,:)./rou;
end
beta_snr_local(ratio_mismatch_num+1,:)=lamda_snrtimes(1,:)*(10.^((SNR(9)-8)/10))/sum(lamda_snrtimes(1,:))./rou;
beta_snr_local(ratio_mismatch_num+2,:)=lamda_snrtimes(1,:)*(10.^((SNR(9)+8)/10))/sum(lamda_snrtimes(1,:))./rou;

r_expand1=zeros(1,site_num);
r_expand2=zeros(1,site_num);
r_expand3=zeros(1,site_num);
r1=zeros(1,trial_num1);
pro_standard=[0.9];
%pro_standard=[0.9 0.7 0.5 0.3];%H1下分为点
expand_num=length(pro_standard);


eta01=zeros(mismatch_num,site_num);
for j=1:mismatch_num
    for i=1:site_num
        randsig=(randn(M(i)+1,trial_num1)+1i*randn(M(i)+1,trial_num1))/sqrt(2);
        loss_factor=betarnd(K(i)-L(i)+2,L(i)-1,1,trial_num1);
        r1=(K(i)+1)*log(1+abs(randsig(1,:)+sqrt(loss_factor)*sqrt(beta_snr_local(j,i))).^2./sum(abs(randsig(2:end,:)).^2));
        r_expand1(i)=mean(r1);
    end
    eta01(j,:)=1-exp(-r_expand1./(K+1));%均值--展开点 取值范围在0~1之间
end
snr_matrixE=beta_snr_local;
file_name=['MismatchMeanValue' '.mat'];
save(file_name,'eta01', 'snr_matrixE')
toc