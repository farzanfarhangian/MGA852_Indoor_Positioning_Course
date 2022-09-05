function zupt = zero_velocity_detector(u,type)
u=u';
sigma_a=0.01;
sigma_g=0.1*pi/180; % 0.1 
Window_size=3; %3
gamma=0.3e5;  %0.3e5
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
%zupt=zupt(end);
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
Window_size=15;
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
sigma_g=1.8*pi/180; %1.6
Window_size=25;
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
