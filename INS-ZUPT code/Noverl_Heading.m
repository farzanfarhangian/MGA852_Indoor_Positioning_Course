function [ori] = Noverl_Heading(u,magno,ori_p)
fb=u(1:3)';  gb=u(4:6)';  mb=magno'; 
Cbn=Rt2b([ori_p(1) ori_p(2) ori_p(3)]);
Cnb=inv(Cbn);     
ggb=Cnb*[0;0;-1];
gbfb = cross(ggb,mb);
gfn = gbfb/norm(gbfb) ; 
mb = cross(gfn,Cnb*[0;0;-1]) ;
q=dcm2q(Cbn);
a=q(1); b=q(2); c=q(3); d=q(4);
H = [c -d a -b;-b -a -d -c;-a b c -d;a b -c -d;-d c b -a;c d a b] ;
q_new=H'*[fb;mb];
R=q2dcm(q_new);

phi=atan2(R(3,2),R(3,3));
theta=-atan(R(3,1)/sqrt(1-R(3,1)^2));
psi=atan2(R(2,1),R(1,1));

ori=[phi;theta;psi];


end