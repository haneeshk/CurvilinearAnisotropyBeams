function [ur, utheta, w]=ComputeDisplacement(kappa_x,kappa_y,theta,r,K,Kp,cV,epsilon,UpVec,UppVec,VpVec,WpVec,mVec,muVec,nu,z)

U=0;
for i=1:4
    U=U+K(i,1)*UpVec(i,1)*r^(mVec(i,1));
end
U=U+UpVec(5,1)*(r^2);
U=(kappa_x*sin(theta)-kappa_y*cos(theta))*U;

for i=1:2
    U=U+Kp(i,1)*UppVec(i,1)*r^(mVec(i,1));
end

U=U+muVec(3,1)*cV*(r^2)+UppVec(4,1)*epsilon*r;

V=0;
for i=1:4
    V=V+K(i,1)*VpVec(i,1)*r^(mVec(i,1));
end
V=V+VpVec(5,1)*(r^2);
V=(kappa_x*cos(theta)+kappa_y*sin(theta))*V;


W=0;
for i=1:4
    W=W+K(i,1)*WpVec(i,1)*r^(mVec(i,1));
end
W=W+WpVec(5,1)*r^2;
W=W*(kappa_x*cos(theta)+kappa_y*sin(theta));


omega1=0;
omega2=0;
omega3=0;
w0p=0;

u0p=-nu*kappa_y;
v0p=nu*kappa_x;

urp=z*(-omega1*sin(theta)+omega2*cos(theta))+u0p*cos(theta)+v0p*sin(theta);
uthetap=z*(-omega1*cos(theta)-omega2*sin(theta))+omega3*r-u0p*sin(theta)+v0p*cos(theta);
wp=r*(omega1*sin(theta)-omega2*cos(theta))+w0p;


ur=-z^2/2*(kappa_x*sin(theta)-kappa_y*cos(theta))+U+urp;

utheta=-z^2*(kappa_x*cos(theta)+kappa_y*sin(theta))+cV*r*z+V+uthetap;

w=z*(kappa_x*r*sin(theta)-kappa_y*r*cos(theta)+epsilon)+W+wp;

end
