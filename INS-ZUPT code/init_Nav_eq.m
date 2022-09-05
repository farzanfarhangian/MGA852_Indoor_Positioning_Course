function [x,q]=init_Nav_eq(u)
f_u=mean(u(1:20,1));
f_v=mean(u(1:20,2));
f_w=mean(u(1:20,3));

roll=atan2(-f_v,-f_w);
pitch=atan2(f_u,sqrt(f_v^2+f_w^2));


% Set the attitude vector
attitude=[roll pitch 150]';

% Calculate quaternion corresponing to the initial attitude
Rb2t=Rt2b(attitude)';
q=dcm2q(Rb2t);

% Set the initial state vector
x=zeros(9,1);
x(1:3,1)=[0;0;0];
x(7:9,1)=attitude;
end

