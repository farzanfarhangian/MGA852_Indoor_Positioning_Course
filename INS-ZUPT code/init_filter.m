function [P,Q,R,H]=init_filter

sigma_initial_acc_bias=0.0005*ones(3,1);        %0.0005*ones(3,1);
sigma_initial_gyro_bias=0.9*pi/180*ones(3,1);   %0.9*pi/180*ones(3,1);
acc_bias_driving_noise=0.0000001*ones(3,1);     % 0.0000001
gyro_bias_driving_noise=0.0000001*pi/180*ones(3,1); 
sigma_initial_pos=1e-2*ones(3,1);            
sigma_initial_vel=1e-2*ones(3,1);               
sigma_initial_att=(pi/180*[0.1 0.1 0]'); 
sigma_acc =0.5*[1 1 1]';
sigma_gyro =0.5*[1 1 1]'*pi/180;
%sigma_vel=[0.0001 0.1];   %real
sigma_vel=[0.0001 0.01 0.01 0.01 0.02];   %compass 0.01 0.0001
%sigma_vel=[0.01 0.01 0.01 0.01 0.01];   %compass altimeter


%  Initial state covariance matrix
P=zeros(9+6);
P(10:12,10:12)=diag(sigma_initial_acc_bias.^2);
P(13:15,13:15)=diag(sigma_initial_gyro_bias.^2);
    
% Process noise covariance matrix
Q=zeros(12);
Q(7:9,7:9)=diag(acc_bias_driving_noise.^2);
Q(10:12,10:12)=diag(gyro_bias_driving_noise.^2);
    
% Observation matrix
%H=zeros(2,9+6); %real
H=zeros(5,9+6);  % compass
%H=zeros(5,9+6);  % compass+altitude

% General values for the observation matrix H
%H(1:3,4:6)=eye(3); %real
%H(3:5,4:6)=eye(3); % compass
%H(2,9)=1; % compass
%H(2,9)=1; % heading
H(1,3)=1;  %altitude
H(2:4,4:6)=eye(3);  %v
H(5,9)=1;  %altitude+ori
%H(6:8,13:15)=eye(3);

% General values for the initial covariance matrix P
P(1:3,1:3)=diag(sigma_initial_pos.^2);
P(4:6,4:6)=diag(sigma_initial_vel.^2);
P(7:9,7:9)=diag(sigma_initial_att.^2);

% General values for the process noise covariance matrix Q
Q(1:3,1:3)=diag(sigma_acc.^2);
Q(4:6,4:6)=diag(sigma_gyro.^2);

% General values for the measurement noise matrix R
R=diag(sigma_vel.^2);

end

