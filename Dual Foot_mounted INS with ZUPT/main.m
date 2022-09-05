%% Loads the algorithm settings and the IMU data
%clear all
disp('Loads the algorithm settings and the IMU data')
[u, heading]=settings(IMU_left', IMU_right', ekf_heading);
u_c = [IMU_right(:,1:3)*9.8 IMU_right(:,4:6)]';

%% Run the zero-velocity detector
global simdata
disp('Runs the zero velocity detector')

% Run the zero-velocity detector for right foot data
[zupt_r T_r]=zero_velocity_detector(u(1:6,:));

% Run the zero-velocity detector for left foot data
[zupt_l T_l]=zero_velocity_detector(u(7:12,:));

% Combine the Right and Left foot zupt data
zupt=[zupt_r ; zupt_l];
clear zupt_*
clear T_*

%% Run the Kalman filter
disp('Runs the filter')
[x_h cov]=ZUPTaidedINS(u, u_c,zupt,heading);

%% Print the horizontal error and spherical error at the end of the
%% trajectory
% sprintf('Horizontal error = %0.5g , Spherical error = %0.5g',sqrt(sum((x_h(1:2,end)).^2)), sqrt(sum((x_h(1:3,end)).^2)))

%sprintf('Magnetic Heading Update = %s , Roll Pitch Update from Accelerometer = %s',simdata.mhupt, simdata.angleupdatefromaccel)

%% View the result 
disp('Views the data')
view_data;

%% Save the workspace
%save([simdata.path 'Workspace.mat'])