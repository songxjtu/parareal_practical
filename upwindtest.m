% Solve the cooled RFR model
% \theta_t = K1 \theta_{xx} -K2 \theta_x +K3 g(\theta)(1- \chi) +K3(1- \theta)
% \chi_t   = K5 \chi_{xx} -K6 \chi_x +K7 g(\theta)(1- \chi)
% where g(\theta) is a nonliear term.
% finite difference method and upwind scheme.


clear;
clc;

K1=6.9393*1.0e-4;
K2=0.1749;
K3=1.5577*1.0e-6;
K4=0.0175;
K5=2.4038*1.0e-3;
K6=174.06;
%K6=K4;
K7=0.01;

N=5;
h=1/(N-1);
T=1;
m=5;
dt=T/m;
e1=ones(N,1);
DIu=spdiags([e1, -2*e1, e1], [-1, 0 ,1], N, N);
DIu(1,1)=DIu(1,1) +1 -K2*h/K1;
DIu(N,N) = DIu(N,N) + 1;

ADu=spdiags([-e1 e1],[-1 0],N,N);
ADu(1,1) = ADu(1,1) - 1 + K2*h/K1;

DIv=spdiags([e1, -2*e1, e1], [-1, 0 ,1], N, N);
DIv(1,1)=DIv(1,1) +1 - 1/K5*K6*h;
DIv(N,N) = DIv(N,N)+1;

ADv=spdiags([-e1 e1],[-1 0],N,N);
ADv(1,1) = ADv(1,1)-1 +  1/K5*K6*h;


Au=K1*DIu/h/h - K2*ADu/h - K4*eye(N,N);
Av=K5*DIv/h/h - K6*ADv/h;

u(1:N,1) = 2.02;
v(1:N,1) = 1;
u(1:N,1) = [1.0001;1.0001;1.0001;1.0001;1];
v(1:N,1) = [0.0300;0.0080;0.0042; 0.0087;0.0286];
u(1:N,1) = [ 1.3236;2.0503;0.9976;0.8601;1.3078];
v(1:N,1) = [ 2.1011;1.9366;0.8006;1.0842;0.0177];

fu=zeros(N,1);

fu(1,1)= (K1/h/h + K2/h)*K2/K1*h;



for i = 1:m
   
    u(1:N,i+1) = (eye(N,N) - dt*Au)\(u(1:N,i)+(K3*g(u(1:N,i)).*(1 - v(1:N,i))...
     +K4 + fu)*dt);
        v(1:N,i+1) = (eye(N,N) - dt*Av)\(v(1:N,i)+(K7*g(u(1:N,i)).*(1 - v(1: ...
    N,i)))*dt);
    
    %u(1:N,i+1) = (eye(N,N) - dt*Au)\(u(1:N,i));
    %v(1:N,i+1) = (eye(N,N) - dt*Av)\(v(1:N,i));

    %     u(1:N,i+1) = (eye(N,N) - dt*Au)\(u(1:N,i)+fu*dt +K4*dt);
    % v(1:N,i+1) = (eye(N,N) - dt*Av)\(v(1:N,i));
end
x=0:h:1;
t=0:dt:T;
subplot(1,2,1)
mesh(t,x,u)

subplot(1,2,2)
mesh(t,x,v)
