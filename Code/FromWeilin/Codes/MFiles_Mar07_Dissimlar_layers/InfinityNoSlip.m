function EI = InfinityNoSlip(b0,bn,angle,flag_mat,flag_sym)

C=ReadCMat(1,angle,flag_mat,flag_sym); % stiffness constants
B=zeros(6,6,1); % reduced stiffness constants
for i=1:6
    for j=1:6
        B(i,j)=C(i,j)-C(i,3)*C(3,j)/C(3,3);
    end
end

a=B(2,2)*B(4,4)-B(2,4)^2;
b=B(2,4)*(2*B(1,4)+B(2,4)+2*B(5,6))-B(4,4)*(B(1,1)+2*B(1,2)+B(2,2)+B(6,6))-B(2,2)*B(5,5)+B(1,4)^2;
c=B(5,5)*(B(1,1)+2*B(1,2)+B(2,2)+B(6,6))-B(5,6)^2;
m1=(-b+sqrt(b^2-4*a*c))/2/a;
m1=sqrt(m1);
m3=(-b-sqrt(b^2-4*a*c))/2/a;
m3=sqrt(m3);
m=[m1 -m1 m3 -m3];
for i=1:4
    g(i)=B(2,4)*m(i)^2+(B(1,4)+B(2,4))*m(i)^1-B(5,6);
    g(i)=g(i)/(B(4,4)*m(i)^2-B(5,5));
    alpha(i)=pi*(C(1,3)+C(2,3)*(m(i)+1)-C(3,4)*g(i)*m(i));
    alpha(i)=alpha(i)/(m(i)+2)/C(3,3);
end
tempMatrix=([-2*B(1,4)-6*B(2,4)+B(5,6), 4*B(4,4)-B(5,5);...
                -B(1,1)-2*B(1,2)+3*B(2,2)-B(6,6), 2*B(1,4)-2*B(2,4)+B(5,6)]);
tempVec=[2*C(3,4), C(1,3)-C(2,3)]'/C(3,3);
mu=tempMatrix\tempVec;
beta=mu(1)*(C(1,3)+3*C(2,3))-2*mu(2)*C(3,4)-1;
beta=pi*beta/C(3,3)/4;

A=[b0^(m(1)-1), b0^(m(2)-1), b0^(m(3)-1), b0^(m(4)-1);...
    g(1)*b0^(m(1)-1), g(2)*b0^(m(2)-1), g(3)*b0^(m(3)-1), g(4)*b0^(m(4)-1);...
    bn^(m(1)-1), bn^(m(2)-1), bn^(m(3)-1), bn^(m(4)-1);...
    g(1)*bn^(m(1)-1), g(2)*bn^(m(2)-1), g(3)*bn^(m(3)-1), g(4)*bn^(m(4)-1)];
y=[-mu(1)*b0, -mu(2)*b0, -mu(1)*bn, -mu(2)*bn]';
KiVec=A\y;

EI=alpha(1)*KiVec(1)*(b0^(2+m(1))-bn^(2+m(1)))...
  +alpha(2)*KiVec(2)*(b0^(2+m(2))-bn^(2+m(2)))...
  +alpha(3)*KiVec(3)*(b0^(2+m(3))-bn^(2+m(3)))...
  +alpha(4)*KiVec(4)*(b0^(2+m(4))-bn^(2+m(4)));

EI=EI+beta*(b0^4-bn^4);

end
