clear;
clc;

K1=6.9393*1.0e-4;
K2=0.1749;
K3=1.5577*1.0e-6;
K4=0.0175;
K5=2.4038*1.0e-3;
K6=174.06;
K7=0.01;

N=40;
h=1/(N-1);
T=1;
m=10;
dt=T/m;
e1=ones(N,1);
DIu=spdiags([e1, -2*e1, e1], [-1, 0 ,1], N, N);
DIu(1,1)=DIu(1,1) - 2/K1*K2*h;
DIu(1,2)=DIu(1,2) +1;
DIu(N,N-1) = DIu(N,N-1)+1;

ADu=spdiags([-e1 e1],[-1 1],N,N);
ADu(1,1) = ADu(1,1) - 2/K1*K2*h;
ADu(1,2)=ADu(1,2) +1;
ADu(N,N-1) = ADu(N,N-1) +1;


DIv=spdiags([e1, -2*e1, e1], [-1, 0 ,1], N, N);
DIv(1,1)=DIv(1,1) - 2/K5*K6*h;
DIv(1,2)=DIv(1,2) +1;
DIv(N,N-1) = DIv(N,N-1)+1;

ADv=spdiags([-e1 e1],[-1 1],N,N);
ADv(1,1) = ADv(1,1) - 2/K5*K6*h;
ADv(1,2)=ADv(1,2) +1;
ADv(N,N-1) = ADv(N,N-1) +1;

Au=K1*DIu/h/h - K2*ADu/2/h + K4*eye(N,N);
Av=K5*DIv/h/h - K6*ADv/2/h;

u(1:N,1) = 2.02;
v(1:N,1) = .7;

fu=zeros(N,1);

fu(1,1)= (K1/h/h + K2/2/h)*K2/K1*2*h;

for i = 1:m
   
   u(1:N,i+1) = (eye(N,N) - dt*Au)\(u(1:N,i)+K3*g(u(1:N,i) + fu).*(1 - v(1:N,i))+K4);
   v(1:N,i+1) = (eye(N,N) - dt*Av)\(v(1:N,i)+K7*g(u(1:N,i)).*(1 - v(1:N,i)));
end

subplot(1,2,1)
mesh(u)

subplot(1,2,2)
mesh(v)
