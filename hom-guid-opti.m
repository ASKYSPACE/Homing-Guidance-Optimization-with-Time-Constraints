clc
clear 
format long;
global Uu V 
P = 7;
V = 380;
r0 = -27500;
t0       = -1;
lambda0  = 0*pi/180;
epsilon0 = 40*pi/180;
rf       = -1;
tf       = -72.95;
lambdaf  = -45*pi/180; 
epsilonf = 5*pi/180;  %If negative, it is easy to have no solution!
R    = r0 - rf;
dr   = 200;
Rr   = r0:dr:rf;
Num    = length(Rr);
sigma0 = tan(epsilon0);
sigmaf = tan(epsilonf);
X0   = [t0;lambda0;sigma0];
Xf   = [tf;lambdaf;sigmaf];
for n = 1:Num
 Uu(n) = 1*interp1([r0,(r0-rf)/2,rf],[3 3 3],Rr(n));      %Part1
end
%% %-------------------------初始序列-------------------%
Uk(:,1) = Uu;
lxk0_1=linspace(X0(1),Xf(1),length(Rr));
lxk0_2=linspace(X0(2),Xf(2),length(Rr));
lxk0_3=linspace(X0(3),Xf(3),length(Rr));
lxk0=[lxk0_1;lxk0_2;lxk0_3];
Xk(:,:,1) = lxk0;
%% %-------------------------边界与约束-------------------%
tmin       = tf-5;            tmax    = t0;
lambdamin  = -80*pi/180;  lambdamax   = 0;
sigmamin   = tan(1*pi/180); sigmamax = tan(80*pi/180); 
lambdafmin = -85*pi/180; lambdafmax   = -45*pi/180;
Nmin      = -6;            Nmax     = 6;
xet = (tmax-tmin)/10;
xelambda = (lambdamax-lambdamin)/10;
xesigma  = (sigmamax-sigmamin)/10; 
ue       = 2.4;
objk(1) = 0;
%%
for s = 1:P-1
       sigma_k = Xk(3,:,s);
    for si = 1:Num
        fk(:,:,si)   =  [-Rr(si)*sqrt(1+sigma_k(si)^2)/V;sigma_k(si);-sigma_k(si)*(1+sigma_k(si)^2)];
        gk(:,:,si)   =  [0;0;sigma_k(si)*(1+sigma_k(si)^2)];        
        Dfk(:,:,si)  =  [0 0 -sigma_k(si)*Rr(si)/(sqrt(1+sigma_k(si)^2)*V);0 0 1;0 0 -1-3*sigma_k(si)^2];    
        Ak(:,:,si)   =  Dfk(:,:,si);
        Bk(:,:,si)   =  gk(:,:,si);
        Ck(:,:,si)   =  fk(:,:,si) - Ak(:,:,si)*Xk(:,si,s);
        Chik(:,:,si) =  Rr(si)^2*sqrt(1 + sigma_k(si)^2)/(V*sigma_k(si)^2);
    end  
    X  = sdpvar(3,Num);
    Aa = sdpvar(1,Num-1);
    U  = sdpvar(1,Num-1);
    F  = [X(:,1) == X0];
    F = [F;X(2,Num) == Xf(2)];
    F = [F;X(3,Num) == Xf(3)]; 
    F = [F; tmin      <= X(1,:) <= tmax];
    F = [F; lambdamin <= X(2,:) <= lambdamax];
    F = [F; sigmamin  <= X(3,:) <= sigmamax];
    F = [F;-xet      <= ((X(1,:)-Xk(1,:,s))) <= xet];
    F = [F;-xelambda <= ((X(2,:)-Xk(2,:,s))) <= xelambda]; 
    F = [F;-xesigma  <= ((X(2,:)-Xk(2,:,s))) <= xesigma];  
    for ii = 1:Num-1
    F =  [F;U(ii)^2   <=  Aa(ii)*Chik(ii)];
    end
    F = [F;Nmin <= U <= Nmax];         
    
    for i = 1:Num-1
        F = [F;X(:,i+1)==dr*(Ak(:,:,i)*X(:,i) + Bk(:,:,i)*U(i) + Ck(:,:,i))/Rr(i)+X(:,i)];
    end
 
    obj = sum(Aa*dr);        
    ops = sdpsettings('solver','ecos','verbose',1);%
    solvesdp(F, obj, ops)
    fprintf('第%d次序列凸优化',s);
    Xktemp = double(X);
    Uk(1:Num-1,s+1) = double(U);
    Uk(Num,s+1) =   Uk(Num-1,s+1);
    Uu = Uk(:,s+1);
    objk(s+1)=double(obj);
    Xk(:,:,s+1) = Xktemp;
    subplot(3,1,1)
    plot(Rr,Uk(:,s+1),'r','LineWidth', 1);xlabel('range of LOS (m)'); ylabel('Gain');
    hold on; plot(Rr,Uk(:,1),'k--','LineWidth', 0.5);
    subplot(3,1,2)
    plot(Rr,Xktemp(1,:),'b'); hold on; 
    subplot(3,1,3)
    plot(Rr,Xktemp(2,:),'b'); hold on; 
end
%%
close all
t = Xktemp(1,:) - Xktemp(1,end);%
for it = 1:P
yr(:,it) = Rr(it).*cos(Xk(2,:,it));
xr(:,it) = Rr(it).*sin(Xk(2,:,it));
end

figure(1)
for j = 1:P
    if j == P
        plot(t,Uk(:,j),'r--','LineWidth', 3);set(gca, 'xdir','reverse');xlabel('time to go (s)'); ylabel('Gain');
    else
       plot(t,Uk(:,j),'LineWidth', 1);set(gca, 'xdir','reverse');xlabel('time to go (s)'); ylabel('Gain');
    end
hold on
end
figure(2)
for j = 1:P
    if j == P
        plot(t,Xk(2,:,j)*180/pi,'r--','LineWidth', 3);set(gca, 'xdir','reverse');xlabel('time to go (s)'); ylabel('\lambda');
    else
        plot(t,Xk(2,:,j)*180/pi,'LineWidth', 1);set(gca, 'xdir','reverse');xlabel('time to go (s)'); ylabel('\lambda');
    end
hold on
end
figure(3)
for j = 1:P
    if j == P
        plot(t,Xk(3,:,j)*180/pi,'r--','LineWidth', 3);set(gca, 'xdir','reverse');xlabel('time to go (s)'); ylabel('\epsilon');
    else
        plot(t,Xk(3,:,j)*180/pi,'LineWidth', 1);set(gca, 'xdir','reverse');xlabel('time to go (s)'); ylabel('\epsilon');
    end
hold on
end
figure(4)
for j = 1:P
    if j == P
        plot(xr(:,j)/1000,-yr(:,j)/1000,'r','LineWidth',3);%set(gca, 'ydir','reverse');xlabel('Range to go (km)'); ylabel('Altitude (m)');
    else
        plot(xr(:,j)/1000,-yr(:,j)/1000,'LineWidth',1); %set(gca, 'ydir','reverse');xlabel('Range to go (km)'); ylabel('Altitude (m)');
    end
hold on
end