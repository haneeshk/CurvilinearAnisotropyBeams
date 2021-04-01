function [U,V,W,ur,utheta,w]=ComputeDisplacement(kappa_x,kappa_y,theta,r,Kn,Kpn,cV,epsilon,Upn,Uppn,Vpn,Wpn,mVecn,muVecn,z,nun)

U=0;
for i=1:4
    U=U+Kn(i,1)*Upn(i,1)*r^(mVecn(i,1));
end
U=U+Upn(5,1)*(r^2);
U=(kappa_x*sin(theta)-kappa_y*cos(theta))*U;

for i=1:2
    U=U+Kpn(i,1)*Uppn(i,1)*r^(mVecn(i,1));
end

U=U+muVecn(3,1)*cV*(r^2)+Uppn(4,1)*epsilon*r;

V=0;
for i=1:4
    V=V+Kn(i,1)*Vpn(i,1)*r^(mVecn(i,1));
end
V=V+Vpn(5,1)*(r^2);
V=(kappa_x*cos(theta)+kappa_y*sin(theta))*V;


W=0;
for i=1:4
    W=W+Kn(i,1)*Wpn(i,1)*r^(mVecn(i,1));
end
W=W+Wpn(5,1)*r^2;
W=W*(kappa_x*cos(theta)+kappa_y*sin(theta));


omega1=0;
omega2=0;
omega3=0;
w0p=0;

u0p=-nun*kappa_y;
v0p=nun*kappa_x;

urp=z*(-omega1*sin(theta)+omega2*cos(theta))+u0p*cos(theta)+v0p*sin(theta);
uthetap=z*(-omega1*cos(theta)-omega2*sin(theta))+omega3*r-u0p*sin(theta)+v0p*cos(theta);
wp=r*(omega1*sin(theta)-omega2*cos(theta))+w0p;




ur=-z^2/2*(kappa_x*sin(theta)-kappa_y*cos(theta))+U+urp;

utheta=-z^2*(kappa_x*cos(theta)+kappa_y*sin(theta))+cV*r*z+V+uthetap;

w=z*(kappa_x*r*sin(theta)-kappa_y*r*cos(theta)+epsilon)+W+wp;






end
