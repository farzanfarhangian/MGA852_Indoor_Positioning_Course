%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------- PDR+EKF indoor positioning algorithm --------------------
%------------- By: Farzan Farhangian -----------------------------------
%------------- January 2021 : LASSENA - ETS University -----------------
% ------------ Copyright ©2022 All rights reserved ---------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data definitions test 1
% Data accel, heading, time ---------------------------------------------
yaw = ekf_heading(1:65000)-110;
accel=IMU_left(1:65000,1:3)*1;
gyro=IMU_left(1:65000,4:6);
time_axis=0:0.01:length(accel)/100;   time_axis=time_axis(1:length(accel));
u=[accel gyro];
q=eul2quat([0 0 (150)*pi/180]); q=q';      %  initial quaternion

%% EKF Parameters
P=diag([1e-6 1e-6]);
Q=diag([1e-7 1e-7]);
R=1.5e-4;

%% Step detection
[StepCount,PkValue,PeakLocation,bandPass_acc] = Stepdetector(accel,99,0.75,2.50,0.14,410,time_axis) ; 

%% PDR/EKF
PositionX = zeros(StepCount, 1);  
PositionY = zeros(StepCount, 1);  
distance = 0;

for k = 1:StepCount-1
    
    % Step length calculator
    [StepLength,distance,heading] = Lengthestimator(accel,yaw,PkValue,PeakLocation,distance,k);
    % position prediction
    PositionX(k+1) = PositionX(k) + StepLength * cos(deg2rad(heading));
    PositionY(k+1) = PositionY(k) + StepLength * sin(deg2rad(heading));
    % EKF Calibration 
    [P,PositionX(k+1),PositionY(k+1)] = EKF_PDR(PositionX(k+1),PositionX(k),PositionY(k+1),PositionY(k),StepLength,P,Q,R);
    
end

%% Error Analysis EKF - UKF
error_ekf=[];
for i=1:length(PositionX)
    d=[];
    pt = [PositionX(i),PositionY(i),0];
    v1=[ref_p(1,2),ref_p(1,1),0];   v2 = [ref_p(2,2),ref_p(2,1),0];
    v3=[ref_p(3,2),ref_p(3,1),0];   v4 = [ref_p(4,2),ref_p(4,1),0];
    
    d(1) = point_to_line(pt, v1, v2);
    d(2) = point_to_line(pt, v2, v3);
    d(3) = point_to_line(pt, v3, v4);
    d(4) = point_to_line(pt, v4, v1);
    
    error_ekf(i)=min(d);
end
e_ekf=mean(error_ekf);

%% Visualization
tiledlayout(2,2) % Structure
title('Detected steps and accel peaks')

% Right plot
nexttile([2 2])
scatter(PositionY,PositionX,'o'); grid on; hold on;% north-up 
plot(ref_p(:,1),ref_p(:,2),'-','Color','#EDB120','LineWidth',1.5);
title('Reference and estimated 2D position')
set(gcf,'position',[50 70 1200 500])
legend('EKF 2D Estimated Position','UKF 2D Estimated Position','Reference')

% Left plot
nexttile([1 2])
plot(bandPass_acc); hold on; plot(PeakLocation, PkValue(:,1),'o');
ylim([-0.8 2])
title('Detected steps and accel peaks')
legend('Bandpass Accel','Steps')
nexttile([1 2])
plot(gyro(:,1)); hold on; plot(gyro(:,2)); hold on; plot(gyro(:,3)); 
xlim([8000 9000])
title('Angular velocity')
legend('Gyto x','Gyro y','Gyro z')

h = msgbox(sprintf('Number of steps: %6g \nTotall distance EKF (m): %2.3g \nTotall time (s): %11g \nEKF Mean Error (m): %6g',...
    StepCount,distance,length(accel(:,1))/100, e_ekf),'Info'); 


%% Attached Functions
function [StepCount,PkValue,PeakLocation,bandPass_acc] = Stepdetector(accel,sampling_freq,freq1,freq2,noise_th,peak_th,time_axis) 
%% preprocessing
acc = sqrt(accel(:,1).^2 + accel(:,2).^2 + accel(:,3).^2);
acc = acc-mean(acc);

%% BANDPASS FILTER
%fs=99;
fs=sampling_freq;
%f1=0.75;               % cuttoff low frequency to get rid of baseline wander
%f2=2.75;                 % cuttoff frequency to discard high frequency noise
%f1=0.75 ; f2=2.75 ;
Wn=[freq1 freq2]/(fs/2);    % cutt off frequency based on fs
N = 4;                % order of 3 less processing
[a,b] = butter(N,Wn); % bandpass filtering
bandPass_acc = filtfilt(a,b,acc);
%bandPass = bandPass/ max(abs(bandPass));

%% find peaks
% Find signal peaks - peaks under a threshold value are considered as noise.

[PkValue, PeakLocation] = findpeaks(bandPass_acc, 'MINPEAKHEIGHT', noise_th); %noise_th = 0.25

%% time interval (steps between 0.4s and 2s)
PkValue(:,2) = time_axis(PeakLocation)*1000;
PkValue(2:end,2) = PkValue(2:end,2)-PkValue(1:end-1,2);
index = find(PkValue(:,2)<peak_th);  %peak_th=400
if isempty(index) == 0
    pos_del = [];
    for k = 2:length(index)
        temp = index(k); % position of the suspicious samples
        if PkValue(temp,1) <= PkValue(temp-1,1)
            pos_del = [pos_del; temp];
        else
            pos_del = [pos_del; temp-1];
        end
    end
    PeakLocation(pos_del) = [];
    PkValue(pos_del,:) = [];
end
StepCount = length(PeakLocation); % step number

end
function [StepLength,distance,heading] = Lengthestimator(accel,yaw,PkValue,PeakLocation,distance_p,k)
    acc = sqrt(accel(:,1).^2 + accel(:,2).^2 + accel(:,3).^2);
    acc = acc-mean(acc);
    pos_start = PeakLocation(k);
    pos_end = PeakLocation(k+1);
    % orientation (yaw)
    %YawSin = mean(yaw(pos_start:pos_end,2));
    %heading=1;
    heading = mean(yaw(pos_start:pos_end));   %correct
    %YawCos = mean(yaw(pos_start:pos_end,3));
    %YawCos = mean(yaw(pos_start:pos_end));
    % step length estimation
    % SL = 0.2844 + 0.2231*frequency + 0.0426*AV
    StepFreq = 1000/PkValue(k+1,2);
    StepAV = var(acc(pos_start:pos_end));
    %StepLength = 0.0844 + 0.6731*StepFreq + 0.14726*StepAV;
    StepLength = 0.2844 + 0.0131*StepFreq + 0.00826*StepAV;       %true:test1
%     StepLength = 0.2844 + 0.131*StepFreq + 0.01026*StepAV;      % test 2
    distance = distance_p + StepLength;
end
%%
function [P,PositionX,PositionY] = EKF_PDR(PositionX_k,PositionX_kp,PositionY_k,PositionY_kp,z,P,Q,R)
F=[1 0;0 1]; % state transition matrix - model the system
x_1=[PositionX_k;PositionY_k];  % state vector
P=F*P*F'+ Q;
P=0.5*(P+P');  %normalizing the P
%H= (1/(sqrt((PositionX_k-PositionX_kp)^2+(PositionY_k-PositionY_kp)^2)))*[x_1(1) x_1(2)];
H= (1/(sqrt((PositionX_k-PositionX_kp)^2+(PositionY_k-PositionY_kp)^2)))*[PositionX_k-PositionX_kp, PositionY_k-PositionY_kp];
K=P*H'*(H*P*H'+R)^(-1);
z_hat=z-H*x_1;
x=x_1+K*(z_hat);
P=(eye(2,2)-K*H)*P;
P=0.5*(P+P');
PositionX=x(1); 
PositionY=x(2);
end
%%
function [ori,q]=Orientation(gyro,q)
Ts=1/100;
w_tb=gyro;
P=w_tb(1)*Ts;
Q=w_tb(2)*Ts;
R=w_tb(3)*Ts;
OMEGA=zeros(4);
OMEGA(1,1:4)=0.5*[0 R -Q P];
OMEGA(2,1:4)=0.5*[-R 0 P Q];
OMEGA(3,1:4)=0.5*[Q -P 0 R];
OMEGA(4,1:4)=0.5*[-P -Q -R 0];
v=norm(w_tb)*Ts;
if v~=0
    q=(cos(v/2)*eye(4)+2/v*sin(v/2)*OMEGA)*q;
    q=q./norm(q);
end
Rb2t=q2dcm(q);
% roll
phi=atan2(Rb2t(3,2),Rb2t(3,3));
%y(7)=ori(1);
% pitch
theta=-atan(Rb2t(3,1)/sqrt(1-Rb2t(3,1)^2));
%y(8)=ori(2);
%yaw
psi=atan2(Rb2t(2,1),Rb2t(1,1));
ori=[phi;theta;psi];
end
function R=q2dcm(q)
p=zeros(6,1);
R=zeros(3,3);
p(1:4)=q.^2;

p(5)=p(2)+p(3);

if p(1)+p(4)+p(5)~=0
   p(6)=2/(p(1)+p(4)+p(5)); 
else
   p(6)=0;
end


R(1,1)=1-p(6)*p(5);
R(2,2)=1-p(6)*(p(1)+p(3));
R(3,3)=1-p(6)*(p(1)+p(2));

p(1)=p(6)*q(1); 
p(2)=p(6)*q(2);
p(5)=p(6)*q(3)*q(4);
p(6)=p(1)*q(2);

R(1,2)=p(6)-p(5);
R(2,1)=p(6)+p(5);

p(5)=p(2)*q(4);
p(6)=p(1)*q(3);

R(1,3)=p(6)+p(5);
R(3,1)=p(6)-p(5);

p(5)=p(1)*q(4);
p(6)=p(2)*q(3);

R(2,3)=p(6)-p(5);
R(3,2)=p(6)+p(5);


end
function q=dcm2q(R)

T = 1 + R(1,1) + R(2,2) + R(3,3);

if T > 10^-8
    
    S = 0.5 / sqrt(T);
    qw = 0.25 / S;
    qx = ( R(3,2) - R(2,3) ) * S;
    qy = ( R(1,3) - R(3,1) ) * S;
    qz = ( R(2,1) - R(1,2) ) * S;

else
    
    if (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))
        
        S = sqrt( 1 + R(1,1) - R(2,2) - R(3,3)) * 2; % S=4*qx
        qw = (R(3,2) - R(2,3)) / S;
        qx = 0.25 * S;
        qy = (R(1,2) + R(2,1)) / S;
        qz = (R(1,3) + R(3,1)) / S;
        
    elseif (R(2,2) > R(3,3))
        
        S = sqrt( 1 + R(2,2) - R(1,1) - R(3,3) ) * 2; %S=4*qy
        qw = (R(1,3) - R(3,1)) / S;
        qx = (R(1,2) + R(2,1)) / S;
        qy = 0.25 * S;
        qz = (R(2,3) + R(3,2)) / S;

    else
        
        S = sqrt( 1 + R(3,3) - R(1,1) - R(2,2) ) * 2; % S=4*qz
        qw = (R(2,1) - R(1,2)) / S;
        qx = (R(1,3) + R(3,1)) / S;
        qy = (R(2,3) + R(3,2)) / S;
        qz = 0.25 * S;

    end

end

%Store in vector
q = [qx qy qz qw]';
end
function zupt = zero_velocity_detector(u,type)
u=u';
sigma_a=0.01;
sigma_g=0.1*pi/180; 
Window_size=3;
gamma=0.3e5;
% Allocate memmory
zupt=zeros(1,length(u));

% Run the desired detector type. Each detector return a vector with their 
% calculated test statistics T. 
switch type
    
    case 'GLRT'
        T=GLRT(u);    
    
    case 'MV'
        T=MV(u);
        
    case 'MAG'
        T=MAG(u);
        
    case 'ARE'
        T=ARE(u);
        
    case 'FARZAN'
        T=FARZAN(u);
        
    otherwise
        disp('The choosen detector type not recognized. The GLRT detector is used')
        T=GLRT(u);
end

% Check if the test statistics T are below the detector threshold. If so, 
% chose the hypothesis that the system has zero velocity 
W=Window_size;
for k=1:length(T)
    if T(k)<gamma
       zupt(k:k+W-1)=ones(1,W); 
    end    
end

% Fix the edges of the detector statistics
T=[max(T)*ones(1,floor(W/2)) T max(T)*ones(1,floor(W/2))];
end
function T=GLRT(u)
sigma_a=0.01;
sigma_g=0.1*pi/180; 
Window_size=5;
gamma=0.3e5;
g=-9.8;
sigma2_a=sigma_a^2;
sigma2_g=sigma_g^2;
W=Window_size;

N=length(u);
T=zeros(1,N-W+1);

for k=1:N-W+1
   
    ya_m=mean(u(1:3,k:k+W-1),2);
    
    for l=k:k+W-1
        tmp=u(1:3,l)-g*ya_m/norm(ya_m);
        T(k)=T(k)+u(4:6,l)'*u(4:6,l)/sigma2_g+tmp'*tmp/sigma2_a;    
    end    
end

T=T./W;

end
function T=MV(u)
sigma_a=0.01;
sigma_g=0.1*pi/180; 
Window_size=15;
gamma=0.3e5;
sigma2_a=sigma_a^2;
W=Window_size;
N=length(u);
T=zeros(1,N-W+1);
for k=1:N-W+1
   
    ya_m=mean(u(1:3,k:k+W-1),2);
    
    for l=k:k+W-1
        tmp=u(1:3,l)-ya_m;
        T(k)=T(k)+tmp'*tmp;    
    end    
end

T=T./(sigma2_a*W);
end
function T=MAG(u)
sigma_a=0.01;
sigma_g=0.1*pi/180; 
Window_size=5;
gamma=0.3e5;
g=-9.8;
sigma2_a=sigma_a^2;
W=Window_size;

N=length(u);
T=zeros(1,N-W+1);

for k=1:N-W+1
    for l=k:k+W-1
        T(k)=T(k)+(norm(u(1:3,l))-g)^2;    
    end    
end

T=T./(sigma2_a*W);
end
function T=ARE(u)
sigma_a=0.01;
sigma_g=1.6*pi/180; 
Window_size=15;
gamma=0.3e5;
sigma2_g=sigma_g^2;
W=Window_size;
N=length(u);
T=zeros(1,N-W+1);

for k=1:N-W+1
    for l=k:k+W-1
        T(k)=T(k)+norm(u(4:6,l))^2;    
    end    
end

T=T./(sigma2_g*W);
end
function T=FARZAN(u)
u=u;
f=u(1:3,:);
g=u(4:6,:);
N=length(u(1,:));
for i=1:N
 c1(i)=norm(f(:,i));
end
for j=1:N
 c3(i)=norm(g(:,i));
end
T=zeros(1,N);
for i=1:N
    if 9<abs(c1(i)) && 11>abs(c1(i))
        if c3(i)<20*pi/180
            T(i)=1;
        end
    end
end


end
function d = point_to_line(pt, v1, v2)
      a = v1 - v2;
      b = pt - v2;
      d = norm(cross(a,b)) / norm(a);
end