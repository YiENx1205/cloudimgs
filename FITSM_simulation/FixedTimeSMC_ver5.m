%QSP, Homogeneous, J.M.wang

clear
close all
clc


% rhoInf=[0.02,0.02,0.02,0.02,0.02,0.02];
% theta=[0.1,0.1,0.1,0.1,0.1,0.1];
% zetaUpper=[0.8,0.8,0.8,0.8,0.8,0.8];
% zetaLower=[0.8,0.8,0.8,0.8,0.8,0.8];

rho0=1;
rhoInf=0.02;
theta=0.2;
zetaUpper=0.8;
zetaLower=0.8;

k1=0.9;
k2=0;
k3=0.02;
k4=0.8;

alpha=5/3;
beta=3/5;
eta=5/3;

q=0.9;
N=6;

% delta=[2,2,2,2,2,2];
delta=2;
h1=[0.5,0.5,0.5,0.5,0.5,0.5];
h2=[0.07,0.07,0.07,0.07,0.07,0.07];
% h=1;
mu1=1;
mu2=0;
% m=[1400,1400,1400,1400,1400,1400];
m=[1,1,1,1,1,1];

L=[2.5,2.5,2.5,2.5,2.5,2.5];
%initial position
px0=39;
% px=[33.22,27.44,21.66,15.88,10.1,4.32];%0��ʼ���
px=[33,27,21,15,9,3];


%initial velocity 
vx0=2;
% vx=[1.8,1.72,1.61,2.795,2.41,1.5];
% vx=[0,0,0,0,0,0];
vx=[2,2,2,2,2,2];


Delta_t=0.01;
t=0;

ux0=0;
ux=[0,0,0,0,0,0];
% intalpha=zeros(6,6000);
for i=1:N
   intalpha(i)=0;
   intbeta(i)=0;
   
end


for j=1:6000
   
    t=t+Delta_t;
    tm(j)=t+Delta_t; 
    %evolution of leader's velocity
    if(j>0&&j<=300)
        vx0=2; 
        ax0=0;
    elseif(j>300&&j<=500)
        vx0=2+2*(j-300)*Delta_t;
        ax0=2;
    elseif(j>500&&j<=1000)
        vx0=6;
        ax0=0;
    elseif(j>1000&&j<=1200)
        vx0=6-2*(j-1000)*Delta_t;
        ax0=-2;
    elseif(j>1200&&j<=1800)
        vx0=2;
        ax0=0;
    elseif(j>1800&&j<=2000)
        vx0=2+(j-1800)*Delta_t;
        ax0=1;
    elseif(j>2000&&j<=2500)
        vx0=4;
        ax0=0;
    elseif(j>2500&&j<=2700)
        vx0=4-(j-2500)*Delta_t;
        ax0=-1;
    else
        vx0=2;
        ax0=0;
    end
    
    



    e(1,j)=mu1*(px0-px(1)-delta-h1(1)*vx(1)-h2(1)*vx(1)^2-L(1))+mu2*(vx(1)-vx0);
    edot(1)=mu1*(vx0-vx(1)-(h1(1)+2*h2(1)*vx(1))/m(1)*ux(1))+mu2/m(1)*(ux(1)-ax0);
    
    for i=2:N
        e(i,j)=mu1*(px(i-1)-px(i)-delta-h1(i)*vx(i)-h2(i)*vx(i)^2-L(i))+mu2*(vx(i)-vx0);
        edot(i)=mu1*(vx(i-1)-vx(i)-(h1(i)+2*h2(i)*vx(i))/m(i)*ux(i))+mu2/m(i)*(ux(i)-ax0);
        
    end
    
    
    for i=1:N
       rho(i,j)=(rho0-rhoInf)*exp(-theta*j*Delta_t)+rhoInf;
       rhodot(i)=theta*(rhoInf-rho0)*exp(-theta*j*Delta_t); 
       a=log(zetaLower*(zetaUpper*rho(i,j)-e(i,j)));
       b=log(zetaUpper*(e(i,j)+zetaLower*rho(i,j)));
       epsilon(i,j)=0.5*log(zetaUpper*(e(i,j)+zetaLower*rho(i,j)))-0.5*log(zetaLower*(zetaUpper*rho(i,j)-e(i,j)));
       C(i)=0.5*1/(e(i,j)+zetaLower*rho(i,j))-0.5*1/(e(i,j)-zetaUpper*rho(i,j));
       rhod(i)=rhodot(i)/rho(i,j);
       epsilondot(i)=C(i)*(edot(i)-e(i,j)*rhodot(i));
    end
    
    
    for i=1:N
        intalpha(i)=intalpha(i)+abs(epsilon(i,j))^alpha*sign(epsilon(i,j))*Delta_t;
        intbeta(i)=intbeta(i)+abs(epsilon(i,j))^beta*sign(epsilon(i,j))*Delta_t;
        s(i)=epsilon(i,j)+k1*intalpha(i)+k2*intbeta(i);
        sdot(i)=epsilondot(i)+k1*abs(epsilon(i,j))^alpha*sign(epsilon(i,j))...
                +k2*abs(epsilon(i,j))^beta*sign(epsilon(i,j));
        
    end
    
    for i=1:N-1
        S(i,j)=q*s(i)-s(i+1);
    end 
        S(N,j)=q*s(N);
    
    for i=1:N-1
       D(i)=-q*C(i)*rhod(i)*e(i,j)+q*k1*abs(epsilon(i,j))^alpha*sign(epsilon(i,j))...
            +q*k2*abs(epsilon(i,j))^beta*sign(epsilon(i,j))-sdot(i+1);
    end
       D(N)=-q*C(N)*rhod(N)*e(N,j)+q*k1*abs(epsilon(N,j))^alpha*sign(epsilon(N,j))...
            +q*k2*abs(epsilon(N,j))^beta*sign(epsilon(N,j));
        
    ux(1)=m(i)*(mu1/((h1(1)+2*h2(1)*vx(1))*mu1-mu2)*(vx0-vx(1))+1/(q*C(1)*((h1(1)+2*h2(1)*vx(1))*mu1-mu2))*D(1)...
              +1/(q*C(1)*((h1(1)+2*h2(1)*vx(1))*mu1-mu2))*(k3*sign(S(1,j))+k4*abs(S(1,j))^eta*sign(S(1,j))));
    for i=2:N
        ux(i)=m(i)*(mu1/((h1(i)+2*h2(i)*vx(i))*mu1-mu2)*(vx(i-1)-vx(i))+1/(q*C(i)*((h1(i)+2*h2(i)*vx(i))*mu1-mu2))*D(i)...
              +1/(q*C(i)*((h1(i)+2*h2(i)*vx(i))*mu1-mu2))*(k3*sign(S(i,j))+k4*abs(S(i,j))^eta*sign(S(i,j))));
    end
     
    
    px0=px0+vx0*Delta_t;
    for i=1:N
       vx(i)=vx(i)+1/m(i)*ux(i)*Delta_t;
%        if vx()
       px(i)=px(i)+vx(i)*Delta_t;
    end
    
    
    p0(j)=px0;
    v0(j)=vx0;
    for i=1:N
       p(i,j)=px(i);
       v(i,j)=vx(i);
       u(i,j)=ux(i);
    end
    
    
    
end


figure(1);
subplot(2,1,1);
plot(tm,p0,'k',tm,p(1,:),tm,p(2,:),tm,p(3,:),tm,p(4,:),'linewidth',1);
axis([0,60,0,220]);
% axis equal;
xlabel('time');
ylabel('position'); 
legend('p0','p1','p2','p3','p4');

% figure(2);
subplot(2,1,2);
% plot(tm,v0,'k',tm,v(1,:),'r',tm,v(2,:),'g',tm,v(3,:),'b',tm,v(4,:),'c','linewidth',1);
plot(tm,v0,'k',tm,v(1,:),tm,v(2,:),tm,v(3,:),tm,v(4,:),'linewidth',1);
axis([0,60,0,7]);
% axis equal;
xlabel('time');
ylabel('speed');
legend('v0','v1','v2','v3','v4');


figure(3);
% subplot(2,1,1);
plot(tm,e(1,:),tm,e(2,:),tm,e(3,:),tm,e(4,:),tm,-zetaLower*rho(1,:),'--k',tm,zetaUpper*rho(1,:),'-.k','linewidth',1);
axis([0,60,-1,1]);
% axis equal;
xlabel('time');
ylabel('error');
legend('e1','e2','e3','e4','-\theta\rho','\theta\rho');

figure(4);
% subplot(2,1,2);
plot(tm,u(1,:),tm,u(2,:),tm,u(3,:),tm,u(4,:),'linewidth',1);
axis([0,60,-3,3]);
% axis equal;
xlabel('time');
ylabel('control input'); 
legend('u1','u2','u3','u4');

% figure(5);
% plot(tm,S(1,:),'r',tm,S(2,:),'g',tm,S(3,:),'b',tm,S(4,:),'c');
% % axis([0,60]);
% % axis equal;
% xlabel('time');
% ylabel('sliding mode'); 
% legend('S1','S2','S3','S4');


% 初始速度图
% figure(6);
% % plot(tm,v0,'k',tm,v(1,:),'r',tm,v(2,:),'g',tm,v(3,:),'b',tm,v(4,:),'c','linewidth',1);
% plot(tm,v0,'k','linewidth',2);
% axis([0,60,0,7]);
% % axis equal;
% xlabel('time');
% ylabel('speed');
% % legend('v0','v1','v2','v3','v4');




