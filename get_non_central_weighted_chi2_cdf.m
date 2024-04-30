function [cdf] = get_non_central_weighted_chi2_cdf(gate,w,k,lambda)
%gate: 1*M x
%w:weights 1*N
%k:freedom vector 1*N
%lambda_s: non-central parameters 1*N
lambda_s=sqrt(lambda);
N=length(w);
gate_num=length(gate);
cdf=zeros(1,gate_num);
for gg=1:gate_num
    phi_func=@(u) 0;
    tmp_func=@(u) 0;
    rho_u=@(u) 1;
    for nn=1:N
        phi_func=@(u) phi_func(u)+k(nn).*atan(w(nn).*u)+lambda_s(nn).^2.*w(nn).*u./(1+w(nn).^2.*u.^2);
        tmp_func=@(u) tmp_func(u)+0.5*(lambda_s(nn)*w(nn).*u).^2./(1+w(nn).^2.*u.^2);
    end
    phi_func=@(u) (phi_func(u)-gate(gg).*u)*0.5;
    
    for nn=1:N
        rho_u=@(u) rho_u(u).*(1+w(nn).^2.*u.^2).^(0.25*k(nn));
    end
    rho_u=@(u) rho_u(u).*exp(tmp_func(u));
    int_func=@(u) sin(phi_func(u))./rho_u(u)./u;
    cdf(gg)=0.5-integral(int_func,0,inf)/pi;
    a=1-exp(-gate/2);%
end
end

