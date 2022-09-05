%% Pedestrian Dead-Reckoning Navigation (PDR)
%  ------------------------------------------------------------------------
%  Authors: Farzan Farhangian, Hamza Benzerrouk plus MATLAB website
%  resources
%  Contact: farzan.farhangian@etsmtl.ca
%  Website: farzanfarhangian.github.io
%  ------------------------------------------------------------------------

%%
NN=length(accel(:,1));
dt=0.01;
omegab=gyro_z;
N(1)=100;EE(1)=100;
d(1)=10; %%% en metre
v(1)=1; %%%% en m/s
phi(1)=0*pi/180;  
Q=zeros(3,3);
Q(1,1)=2.5;Q(2,2)=2.5;Q(3,3)=(1.5*pi/180)^2;

R=zeros(3,3);
R(1,1)=0.05;R(2,2)=0.05;R(3,3)=(0.5*pi/180)^2;Rmap=R(3,3);

P=zeros(3,3);
P(1,1)=100;P(2,2)=100;P(3,3)=5*pi/180;
P11(1)=P(1,1);P22(1)=P(2,2);P33(1)=P(3,3);

K=zeros(3,3);
H=eye(3,3);
K=P*H'*(H*P*H'+R)^(-1);
K11(1)=K(1,1);K22(1)=K(2,2);%K33(1)=K(3,3);

phiest(1)=phi(1);

Nest(1)=N(1);Eest(1)=EE(1);

W=eye(3,3);

Xest(:,1)=0.95*[Nest(1);Eest(1);phiest(1)];

R0=R;

V(:,1)=[sqrt(R(1,1))*randn;sqrt(R(2,2))*randn;sqrt(R(3,3))*randn]; %%%% observation noise 
VV(:,1)=V(:,1)'*V(:,1);
INOV=H*P*H'+R;
dV=(VV(:,1)*inv(INOV));
alpha=1;
c=3;
  
RR=R;
  

for k=2:NN
    v(k)=v(k-1);
    vb(k)=v(k)+0.01*randn;
%% transition matrix    
    F=[1 0 -dt*vb(k-1)*sin(Xest(3,k-1));...
       0 1 dt*vb(k-1)*cos(Xest(3,k-1));...
       0 0 1];
%% prediction step 
    phiest(k)=Xest(3,k-1)+omegab(k-1)*dt;
    if phiest(k)>=2*pi 
        phiest(k)=phi(k)-2*pi;
    end  
    
    Nest(k)=Xest(1,k-1)+dd(k-1)*cos(Xest(3,k-1));
    Eest(k)=Xest(2,k-1)+dd(k-1)*sin(Xest(3,k-1));

    Xest(:,k-1)=[Nest(k);Eest(k);phiest(k)];

%% EKF update and measurement model
    P=F*P*F'+W*Q*W';
    P=0.5*(P+P');

    % Measurement model
    H=eye(3,3);
    z(1,k)=N(k)+sqrt(R(1,1))*randn;
    z(2,k)=EE(k)+sqrt(R(2,2))*randn;
    z(3,k)=phi(k)+sqrt(R(3,3))*randn;
    z(:,k)=[z(1,k);z(2,k);z(3,k)];

    % Kalman gain
    K=P*H'*(H*P*H'+R)^(-1);
    K11(k)=K(1,1);
    K22(k)=K(2,2);
    K33(k)=K(3,3);

    % Update
    Xest(:,k)=Xest(:,k-1)+K*(z(:,k)-H*Xest(:,k-1));
    Xest1(k)=Xest(1,k);
	Xest2(k)=Xest(2,k);
    Xest3(k)=Xest(3,k);
 
    if Xest3(k)>=2*pi 
        Xest3(k)=xest3(k)-2*pi;
    end
 
    P=(eye(3,3)-K*H)*P;
    P=0.5*(P+P');
    P11(k)=P(1,1);
    P22(k)=P(2,2);
    P33(k)=P(3,3);
end