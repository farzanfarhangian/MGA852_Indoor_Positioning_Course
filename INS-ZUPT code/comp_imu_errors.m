function [u_out]=comp_imu_errors(u_in,x_h)
u_out = u_in + x_h(10:end); %real
%u_out_1=u_in(1:6)+x_h(10:end);  %compass
%u_out = [u_out_1  u_in(7:9)];    % compass
end

