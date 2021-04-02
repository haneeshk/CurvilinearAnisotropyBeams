function [sigmaz, sigmar, sigmatheta, taurtheta, taurz, tauthetaz]=ComputeStressNew(r,theta,kappa_x,kappa_y,epsilon,cV,beta,C,mVec,mpVec,muVec,gVec,gpVec,K,Kp)

sigmar=(kappa_x*sin(theta)-kappa_y*cos(theta))*(K(1)*power(r,mVec(1)-1)+K(2)*power(r,mVec(2)-1)+K(3)*power(r,mVec(3)-1)+K(4)*power(r,mVec(4)-1)+muVec(1)*r)...
            +Kp(1)*power(r,mpVec(1)-1)+Kp(2)*power(r,mpVec(2)-1)+muVec(3)*cV*r+muVec(5)*epsilon;

sigmatheta=(kappa_x*sin(theta)-kappa_y*cos(theta))*(K(1)*(mVec(1)+1)*power(r,mVec(1)-1)+K(2)*(mVec(2)+1)*power(r,mVec(2)-1)+K(3)*(mVec(3)+1)*power(r,mVec(3)-1)+K(4)*(mVec(4)+1)*power(r,mVec(4)-1)+3*muVec(1)*r)...
            +Kp(1)*mpVec(1)*power(r,mpVec(1)-1)+Kp(2)*mpVec(2)*power(r,mpVec(2)-1)+2*muVec(3)*cV*r+muVec(5)*epsilon;
  
taurtheta=-(kappa_x*cos(theta)+kappa_y*sin(theta))*((K(1)*power(r,mVec(1)-1)+K(2)*power(r,mVec(2)-1)+K(3)*power(r,mVec(3)-1)+K(4)*power(r,mVec(4)-1)+muVec(1)*r));

taurz=(kappa_x*cos(theta)+kappa_y*sin(theta))*((K(1)*gVec(1)*power(r,mVec(1)-1)+K(2)*gVec(2)*power(r,mVec(2)-1)+K(3)*gVec(3)*power(r,mVec(3)-1)+K(4)*gVec(4)*power(r,mVec(4)-1)+muVec(2)*r));

tauthetaz=-(kappa_x*sin(theta)-kappa_y*cos(theta))*((K(1)*gVec(1)*mVec(1)*power(r,mVec(1)-1)+K(2)*gVec(2)*mVec(2)*power(r,mVec(2)-1)+K(3)*gVec(3)*mVec(3)*power(r,mVec(3)-1)+K(4)*gVec(4)*mVec(4)*power(r,mVec(4)-1)+2*muVec(2)*r))...
            -(Kp(1)*gpVec(1)*power(r,mpVec(1)-1)+Kp(2)*gpVec(2)*power(r,mpVec(2)-1))-muVec(4)*cV*r-(muVec(5)*(beta(1,4)+beta(2,4))/beta(4,4)+C(3,4)/(C(3,3)*beta(4,4)))*epsilon;

sigmaz=(kappa_x*r*sin(theta)-kappa_y*r*cos(theta)+epsilon-C(1,3)*sigmar-C(2,3)*sigmatheta-C(3,4)*tauthetaz)/C(3,3);

end
