function [F,G]=state_matrix(q,u)
u=u';
Ts = 0.01;
Rb2t=q2dcm(q);
f_t=Rb2t*u(1:3);

% Create a ske symmetric matrix of the specific fore vector
St=[0 -f_t(3) f_t(2); f_t(3) 0 -f_t(1); -f_t(2) f_t(1) 0];

O=zeros(3);
I=eye(3);

acc_bias_instability_time_constant_filter=inf;
gyro_bias_instability_time_constant_filter=inf;
% Correlation constant for accelerometer and gyro biases
B1=-1/acc_bias_instability_time_constant_filter*eye(3);
B2=-1/gyro_bias_instability_time_constant_filter*eye(3);

% Fc=[O I O O O;
%     O O St Rb2t O;
%     O O O O -Rb2t;
%     O O O B1 O;
%     O O O O B2];
Fc=[O I O O O;
    O O -St Rb2t O;
    O O O O Rb2t;
    O O O B1 O;
    O O O O B2];

    
% Noise gain matrix
Gc=[O O O O; Rb2t O O O; O -Rb2t O O; O O I O; O O O I];

F=eye(size(Fc))+Ts*Fc;
G=Ts*Gc;
end

