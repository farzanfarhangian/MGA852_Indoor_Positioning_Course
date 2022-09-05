%% IMU setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%u=[accel(:,1)*9.8 accel(:,2)*9.8 accel(:,3)*9.8 gyro(:,1)*pi/180 gyro(:,2)*pi/180 gyro(:,3)*pi/180];
u=[accel(:,1) accel(:,2) accel(:,3) gyro(:,1) gyro(:,2) gyro(:,3)];
u_f=[accel_f(:,1) accel_f(:,2) accel_f(:,3) gyro_f(:,1) gyro_f(:,2) gyro_f(:,3)];
%u=[accel(:,1)*9.8 accel(:,2)*9.8 accel(:,3)*9.8 gyro(:,1)*pi/180 gyro(:,2)*pi/180 gyro(:,3)*pi/180 magno_calib(:,1) magno_calib(:,2) magno_calib(:,3)];
%u=[accel(:,1) accel(:,2) accel(:,3) gyro(:,1)*pi/180 gyro(:,2)*pi/180 gyro(:,3)*pi/180 magno(:,1) magno(:,2) magno(:,3)];
N=length(u(:,1));

%% ZUPT setting and types %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%zupt=zero_velocity_detector(u_f,'ARE'); %GLRT   MV   MAG   ARE  FAR

%% Initialization for navigation and states %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zupt=[];
[P,Q,R,H]=init_filter; %done
x_h=zeros(N,15);  %[dr,dv,domega,ba,bg]
zupt_vector = zeros(N,1);
[x_h(1:9,1),q]=init_Nav_eq(u); %done
%% INS / ZUPT Equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zupt(1)=0;
for k=2:N 
    
    if k<16
        zupt(k)=0;
    else
        zupt=zero_velocity_detector(u_f(k-15:k,:),'ARE');
    end
    [u_h]=comp_imu_errors(u(k,:),x_h(k-1,:));  %done
    %[x_h(k,:),q]=Navigation_equations(new_ori_pi(k,:)*pi/180,x_h(k-1,:),u_h,q);  %done
    [x_h(k,:),q]=Navigation_equations(x_h(k-1,:),u_h,q);
    [F,G]=state_matrix(q,u_h);  %done
    P=F*P*F'+G*Q*G';
    P=(P+P')/2;
    x_h(k,8)=0;
    
    if zupt(k)==true 
        x_h(k,8)=0;
        %[x_h(k,7:9)] = Noverl_Heading(u(k,:),magno(k,1:3),x_h(k,7:9));
        K=(P*H')/(H*P*H'+R);
        mk=-0.54881;%HDR(x_h(k,9),x_h(k-1,9));
        mk=0;
        %compass = compute_heading_compass(u(k,:),x_h(k-1,:),magno(k,1:3));
        %z=-x_h(k,4:6);  %real
        %z=[-x_h(k,9)+compass -x_h(k,4:6)];  %compass
        %z=[-x_h(k,3) -x_h(k,9)+compass -x_h(k,4:6)];  %compass & true altimeter
        %z=[-x_h(k,4:6) heading(k)-x_h(k,9)];    %heading
        %z=[-x_h(k,3) -x_h(k,4:6) mk x_h(k,9)-(ref(k)*pi/180)];    %mk x_h(k,9)-(ref(k)*pi/180)
        z=[-x_h(k,3) -x_h(k,4:6) 0];  %x_h(k,9)-reff_h(k)
        %z=[-x_h(k,3) -x_h(k,4:6) accel(k,1:3) magno(k,1:3)];
        %z=[-x_h(k,3) -x_h(k,4:6)];
        dx=K*z';
        [x_h(k,:),q]=comp_internal_states(x_h(k,:),dx',q); %done
        P=(eye(15)-K*H)*P;
        P=(P+P')/2;
        zupt_vector(k)=1;
        
        % heading correction new :: TEST %%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            if k>=800
%              x=sort(find(zupt_vector==1),'descend');
%              y=sort(find(zupt_vector==1),'descend');
%              zupt_k=x(1);  zupt_k_1=x(2); zupt_k_2=x(30);
%               if abs(x_h(k,9)-x_h(zupt_k_1,9))<=80*pi/180 &&  abs(x_h(k,9)-x_h(zupt_k_1,9))>=0*pi/180
%                x_h(k,9) = x_h(zupt_k_1,9);
% %                x_h(k,9) = x_h(y(200),9);
% %                elseif abs(x_h(k,9)-x_h(k-2,9))>77*pi/180
% %                   if x_h(k,9)-x_h(k-2,9)<=0
% %                       x_h(k,9)=x_h(k-2,9)-90*pi/180;
% %                   else
% %                       x_h(k,9)=x_h(k-2,9)+90*pi/180;
% %                   end
%               end
%            end
%         % -----------------------------------------------------
% %     else
% %         if k>=800
% %              y=sort(find(zupt_vector==0),'descend');
% %              zupt_k=y(1);  zupt_k_1=y(20); zupt_k_2=y(30); 
% %               if abs(x_h(k,9)-x_h(k-5,9))<=80*pi/180
% %                x_h(k,9) = x_h(k-5,9);
% %               end
% %         end
%      end
%     
%                if k>=800
%               x=sort(find(zupt_vector==1),'descend');
%               y=sort(find(zupt_vector==1),'descend');
%               zupt_k=y(1);  zupt_k_1=y(2); zupt_k_2=y(30);
%                if abs(x_h(k,9)-x_h(zupt_k_1,9))<=20*pi/180
%                    x_h(k,9)=x_h(zupt_k_1,9);
%               end
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
 plot(x_h(:,2),x_h(:,1)); hold on;
 %plot(ref_p(:,1),ref_p(:,2));
 %plot3(x_h(:,1),x_h(:,2),x_h(:,3)); hold on;
 
% show the z-axis down
%  figure(4);
%  plot(time_axis,x_h(:,3)); hold on;
 figure(347);
%plot(time_axis,zupt_vector*360-180); hold on;
plot(time_axis,x_h(:,9)*180/pi);  hold on;
% show estimated velocity
  figure(5);
  plot(time_axis,x_h(:,4)); hold on;
  plot(time_axis,x_h(:,5)); hold on;

% show estimated omega
%  figure(99);
%   plot(time_axis,x_h(:,7)*180/pi); hold on;
%  figure(100);
%   plot(time_axis,x_h(:,8)*180/pi); hold on;
% 
%  figure(101);
%   plot(time_axis,x_h(:,10)); hold on;
%   plot(time_axis,x_h(:,11)); hold on;
% 
%   figure(102);
%   plot(time_axis,x_h(:,13)); hold on;
%   plot(time_axis,x_h(:,14)); hold on;
