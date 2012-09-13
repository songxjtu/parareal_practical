clear
clc

imax=10;
imin=2;
i=imin:imax;
%m=[10 20 40 60 160 320 640];
m=2.^i;
n=2.^i+1;
num=length(n);

[uend, vend] = upwindfunc(m(num), n(num));
erru=[];
errv=[];
for i = 1:num-1  
    [u v] = upwindfunc(m(i), n(i));
    
    erru=[erru max(max(abs(u - uend(1:2^(num-i):end,1:2^(num - i):end))))]
    errv=[errv max(max(abs(v - vend(1:2^(num-i):end,1:2^(num ...
                                                     - i):end))))]
    pause;
end
subplot(1,2,1)
loglog(m(1:end-1),erru,'-',m(1:end-1),1./m(1:end-1),'--')

subplot(1,2,2)
loglog(m(1:end-1),errv,'-',m(1:end-1),1./m(1:end-1),'--')


