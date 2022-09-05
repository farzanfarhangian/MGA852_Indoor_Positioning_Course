%% IMU setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_f=[left_IMU(:,1)*1 left_IMU(:,2)*1 left_IMU(:,3)*1 left_IMU(:,4) left_IMU(:,5) left_IMU(:,6)];
%u_f=[IMU_right(:,1)*9.8 IMU_right(:,2)*9.8 IMU_right(:,3)*9.8 IMU_right(:,4) IMU_right(:,5) IMU_right(:,6)];
%u=[accel(:,1) accel(:,2) accel(:,3) gyro(:,1) gyro(:,2) gyro(:,3)];
u=[chest_IMU(:,1)*1 chest_IMU(:,2)*1 chest_IMU(:,3)*1 chest_IMU(:,4) chest_IMU(:,5) chest_IMU(:,6)];
%u=u_f;
%u=[accel(:,1)*9.8 accel(:,2)*9.8 accel(:,3)*9.8 gyro(:,1)*pi/180 gyro(:,2)*pi/180 gyro(:,3)*pi/180 magno_calib(:,1) magno_calib(:,2) magno_calib(:,3)];
%u=[accel(:,1) accel(:,2) accel(:,3) gyro(:,1)*pi/180 gyro(:,2)*pi/180 gyro(:,3)*pi/180 magno(:,1) magno(:,2) magno(:,3)];
N=length(u(:,1));
heading=ekf_heading(1000:N);
u_f=u_f(1000:N,:);
u=u(1000:N,:);
N=length(u(:,1));

%% ZUPT setting and types %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zupt=zero_velocity_detector(u_f,'ARE'); %GLRT   MV   MAG   ARE  FAR

%% Initialization for navigation and states %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[P,Q,R,H]=init_filter; %done
x_h=zeros(N,15);  %[dr,dv,domega,ba,bg]
zupt_vector = zeros(N,1);
[x_h(1:9,1),q]=init_Nav_eq(u); %done
x_h(9,1)=150*pi/180;
%% INS / ZUPT Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=2:N 
    
    [u_h]=comp_imu_errors(u(k,:),x_h(k-1,:));  %done
    %[x_h(k,:),q]=Navigation_equations(new_ori_pi(k,:)*pi/180,x_h(k-1,:),u_h,q);  %done
    [x_h(k,:),q]=Navigation_equations(x_h(k-1,:),u_h,q);
    [F,G]=state_matrix(q,u_h);  %done
    P=F*P*F'+G*Q*G';
    P=(P+P')/2;
    %x_h(k,7)=0;
    
    if zupt(k)==true 
        %x_h(k,7)=0;
        %[x_h(k,7:9)] = Noverl_Heading(u(k,:),magno(k,1:3),x_h(k,7:9));
        %hhh=Noverl_Heading(u(k,:),magno(k,1:3),x_h(k,7:9));
        K=(P*H')/(H*P*H'+R);
        mk=-0.54881;%HDR(x_h(k,9),x_h(k-1,9));
        mk=0;
        %compass = compute_heading_compass(u(k,:),x_h(k-1,:),magno(k,1:3));
        %z=-x_h(k,4:6);  %real
        %z=[-x_h(k,9)+compass -x_h(k,4:6)];  %compass
        %z=[-x_h(k,3) x_h(k,9)-heading(k)*pi/180];  %compass & true altimeter
        %z=[-x_h(k,4:6) heading(k)-x_h(k,9)];    %heading
        %z=[-x_h(k,3) -x_h(k,4:6) mk x_h(k,9)-(ref(k)*pi/180)];    %mk x_h(k,9)-(ref(k)*pi/180)
        z=[-x_h(k,3) -x_h(k,4:6) x_h(k,9)-(heading(k)+0)*pi/180];  %x_h(k,9)-heading(k)*pi/180
        %z=[-x_h(k,3) -x_h(k,4:6) accel(k,1:3) magno(k,1:3)];
        %z=[-x_h(k,3) -x_h(k,4:6)];
        dx=K*z';
        [x_h(k,:),q]=comp_internal_states(x_h(k,:),dx',q); %done
        P=(eye(15)-K*H)*P;
        %P=(P+P')/2;
        zupt_vector(k)=1;  
    end
end

%% Showing the results and plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_axis=zeros(N,1);
for i=1:N
    time_axis(i,1)=(i-1)/100;
end

% show zupt detecting
figure(1);
plot(time_axis,zupt_vector); hold on; ylim([-0.1 1.1]);

% show the attitude
 figure(2);
  %plot(time_axis,x_h(:,7)*180/pi); hold on;
  %plot(time_axis,x_h(:,8)*180/pi); hold on;
  plot(time_axis,x_h(:,9)*180/pi); hold on;
  %plot(time_axis,ref);
%  plot(time_axis,new_ori(:,3)+150); hold on;
%plot(time_axis,orientation(:,3)); hold on;  %compass
 
 % show the x-y position
 figure(3);
 plot(x_h(10000:N,2)-x_h(10000,2),x_h(10000:N,1)-x_h(10000,1),'LineWidth',1); hold on;
 plot(ref_p(:,1),ref_p(:,2),'r');
 %plot3(x_h(:,1),x_h(:,2),x_h(:,3)); hold on;
 
% show the z-axis down
%  figure(4);
%  plot(time_axis,x_h(:,3)); hold on;
 figure(347);
%plot(time_axis,zupt_vector*360-180); hold on;
plot(time_axis,x_h(:,9)*180/pi);  hold on;
% plot(time_axis,x_h(:,8)*180/pi);  hold on;
% plot(time_axis,x_h(:,7)*180/pi);  hold on;
% show estimated velocity
  figure(5);
  plot(time_axis,x_h(:,4)); hold on;
  plot(time_axis,x_h(:,5)); hold on;

% show estimated omega
%  figure(99);
%   plot(time_axis,x_h(:,7)*180/pi); hold on;
%  figure(100);
%   plot(time_axis,x_h(:,8)*180/pi); hold on;

figure(101);
plot(time_axis,x_h(:,10)); hold on;
plot(time_axis,x_h(:,11)); hold on;
plot(time_axis,x_h(:,12));
% 
figure(102);
plot(time_axis,x_h(:,13)); hold on;
plot(time_axis,x_h(:,14)); hold on;
plot(time_axis,x_h(:,15));

error_ekf=[];
for i=1:length(x_h(10000:N,1))
    d=[];
    pt = [x_h(i,1)-x_h(10000,1),x_h(i,2)-x_h(10000,2),0];
    v1=[ref_p(1,2),ref_p(1,1),0];   v2 = [ref_p(2,2),ref_p(2,1),0];
    v3=[ref_p(3,2),ref_p(3,1),0];   v4 = [ref_p(4,2),ref_p(4,1),0];
    
    d(1) = point_to_line(pt, v1, v2);
    d(2) = point_to_line(pt, v2, v3);
    d(3) = point_to_line(pt, v3, v4);
    d(4) = point_to_line(pt, v4, v1);
    
    error_ekf(i)=min(d);
end
e_ekf=mean(error_ekf);