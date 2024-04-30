close all;clear;
load('PFA4S1762.mat')
pd_expand41=pd_expand_r1;
load('PFA4S8421.mat')
pd_expand42=pd_expand_r1;
p3delta=0.9974;
sig_1=qfuncinv(p3delta/2);
figure;
semilogy(sigma2,pd_expand41,'-+','linewidth',1.5,'markersize',8);hold on;
semilogy(sigma2,pd_expand42,'-o','linewidth',1.5,'markersize',8);hold on;
% semilogy(sigma2,(10^(-4)+sig_1)*ones(1,length(sigma2)),'-','linewidth',1.5,'markersize',8); hold on;
% semilogy(sigma2,(10^(-4)-sig_1)*ones(1,length(sigma2)),'-.','linewidth',1.5,'markersize',8); 

xlim([1 10])
ylim([10^(-5) 10^(-3)])
grid on;
xlabel('\sigma^2');ylabel('Probability of false alarm')
legend('Parameter 1','Parameter 2');
set(gcf,'color',[1,1,1]);
set(gca,'Fontname','Times New Roman','FontSize',13);

load('PFA3S1762.mat')
pd_expand41=pd_expand_r1;
load('PFA3S8421.mat')
pd_expand42=pd_expand_r1;
figure;
semilogy(sigma2,pd_expand41,'-+','linewidth',1.5,'markersize',8);hold on;
semilogy(sigma2,pd_expand42,'-o','linewidth',1.5,'markersize',8);hold on;
xlim([1 10])
ylim([10^(-4) 10^(-2)])
grid on;
xlabel('\sigma^2');ylabel('Probability of false alarm')
legend('Parameter 1','Parameter 2');
set(gcf,'color',[1,1,1]);
set(gca,'Fontname','Times New Roman','FontSize',13);
