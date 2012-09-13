function y = periodicfun(x)
    
    n=length(x)/2;
    u0=x(1:n);
    v0=x(n+1:end);
    m=n;
    [u,v] = upwindfun(u0,v0,m,n);
    y=norm([u0;v0] - [u(end:-1:1,end);v(end:-1:1,end)])
    

