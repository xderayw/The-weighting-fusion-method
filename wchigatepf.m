function x0=wchigatepf(x0,v,w,nint,pf,step)
%根据起始点x0搜索加权卡方分布的接近pf的门限并返回
%   x0=wchigatepf(x0,v,w,nint,pf,step)
%

xpf=wchicdf(x0,v,w,nint)-pf;
if xpf>0
    while xpf>0
        x0=x0+step;
        xpf=wchicdf(x0,v,w,nint)-pf;
    end
elseif xpf<0
    while xpf<0
        x0=x0-step;
        xpf=wchicdf(x0,v,w,nint)-pf;
    end
else
    gate=x0;
    err=0;
    return;
end