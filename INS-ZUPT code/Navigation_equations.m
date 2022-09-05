function [y,q]=Navigation_equations(x,u,q)
u=u';
x=x';
y=zeros(size(x));
Ts=1/99.3393;

%*************************************************************************%
% Update the quaternion vector "q"  given the angular rate measurements.
%*************************************************************************%

w_tb=u(4:6);

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

% Get the roll, pitch and yaw
Rb2t=q2dcm(q);
% roll
y(7)=atan2(Rb2t(3,2),Rb2t(3,3));
%y(7)=ori(1);
% pitch
y(8)=-atan(Rb2t(3,1)/sqrt(1-Rb2t(3,1)^2));
%y(8)=ori(2);
%yaw
y(9)=atan2(Rb2t(2,1),Rb2t(1,1));
%y(9)=ori(3);

%*************************************************************************%
% Update position and velocity states using the measured specific force,
% and the newly calculated attitude.
%*************************************************************************%

% Gravity vector
g_t=[0 0 -9.8]';

% Transform the specificforce vector into navigation coordinate frame.
f_t=q2dcm(q)*u(1:3);

% Subtract (add) the gravity, to obtain accelerations in navigation
% coordinat system.
acc_t=f_t+g_t;

% State space model matrices
A=eye(6);
A(1,4)=Ts;
A(2,5)=Ts;
A(3,6)=Ts;

%B=[(Ts^2)/2*eye(3);Ts*eye(3)];
B=[a*0.5*Ts*Ts*eye(3);b*Ts*eye(3)];    %5.76*[0.5*Ts*Ts*eye(3);Ts*eye(3)];

% Update the position and velocity estimates.
y(1:6)=A*x(1:6)+B*acc_t;
end

