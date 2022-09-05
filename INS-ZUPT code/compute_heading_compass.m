function compass = compute_heading_compass(u,x_h,magno)
u=u';
x_h=x_h';
m_b=magno';
phi=x_h(7); theta=x_h(8);
Hx=m_b(1)*cos(theta)+m_b(2)*sin(phi)*sin(theta)+m_b(3)*cos(phi)*sin(theta);
Hy=m_b(2)*cos(phi)-m_b(3)*sin(phi);
M_d = 120*pi/180;
compass=atan(Hy/Hx)-M_d;
end

