function [sigmaz, sigmar, sigmatheta, taurtheta, taurz, tauthetaz]=ComputeStress(r,theta,kappa_x,kappa_y,epsilon,beta, C, mVec, muVec, gVec, K)

sigmaz=(epsilon - r*kappa_y*cos(theta) + r*kappa_x*sin(theta) - C(2,3)*(-(kappa_y*cos(theta)) + kappa_x*sin(theta))*(K(6,1)/r + power(r,-1 + mVec(1,1))*K(1,1)*(1 + mVec(1,1)) + power(r,-1 + mVec(2,1))*K(2,1)*(1 + mVec(2,1)) + power(r,-1 + mVec(3,1))*K(3,1)*(1 + mVec(3,1)) + power(r,-1 + mVec(4,1))*K(4,1)*(1 + mVec(4,1)) + 3*r*muVec(1,1)) - ....
     C(1,3)*(((-(kappa_y*cos(theta)) + kappa_x*sin(theta))*(K(5,1) + K(6,1) + K(6,1)*log(r) + (power(r,mVec(1,1))*K(1,1)*(1 + mVec(1,1)))/mVec(1,1) + (power(r,mVec(2,1))*K(2,1)*(1 + mVec(2,1)))/mVec(2,1) + (power(r,mVec(3,1))*K(3,1)*(1 + mVec(3,1)))/mVec(3,1) + (power(r,mVec(4,1))*K(4,1)*(1 + mVec(4,1)))/mVec(4,1) + (3*power(r,2)*muVec(1,1))/2.))/r +.... 
        ((kappa_y*cos(theta) - kappa_x*sin(theta))*(r*K(5,1) + r*K(6,1)*log(r) + (power(r,1 + mVec(1,1))*K(1,1))/mVec(1,1) + (power(r,1 + mVec(2,1))*K(2,1))/mVec(2,1) + (power(r,1 + mVec(3,1))*K(3,1))/mVec(3,1) + (power(r,1 + mVec(4,1))*K(4,1))/mVec(4,1) + (power(r,3)*muVec(1,1))/2.))/power(r,2)) + ....
     C(3,4)*(-(kappa_y*cos(theta)) + kappa_x*sin(theta))*(power(r,-1 + mVec(1,1))*gVec(1,1)*K(1,1)*mVec(1,1) + power(r,-1 + mVec(2,1))*gVec(2,1)*K(2,1)*mVec(2,1) + power(r,-1 + mVec(3,1))*gVec(3,1)*K(3,1)*mVec(3,1) + power(r,-1 + mVec(4,1))*gVec(4,1)*K(4,1)*mVec(4,1) + 2*r*muVec(2,1)))/C(3,3);




 
sigmar= ((-(kappa_y*cos(theta)) + kappa_x*sin(theta))*(K(5,1) + K(6,1) + K(6,1)*log(r) + (power(r,mVec(1,1))*K(1,1)*(1 + mVec(1,1)))/mVec(1,1) + (power(r,mVec(2,1))*K(2,1)*(1 + mVec(2,1)))/mVec(2,1) + (power(r,mVec(3,1))*K(3,1)*(1 + mVec(3,1)))/mVec(3,1) + ....
        (power(r,mVec(4,1))*K(4,1)*(1 + mVec(4,1)))/mVec(4,1) + (3*power(r,2)*muVec(1,1))/2.))/r + ((kappa_y*cos(theta) - kappa_x*sin(theta))*....
      (r*K(5,1) + r*K(6,1)*log(r) + (power(r,1 + mVec(1,1))*K(1,1))/mVec(1,1) + (power(r,1 + mVec(2,1))*K(2,1))/mVec(2,1) + (power(r,1 + mVec(3,1))*K(3,1))/mVec(3,1) + (power(r,1 + mVec(4,1))*K(4,1))/mVec(4,1) + (power(r,3)*muVec(1,1))/2.))/power(r,2);
 


      
      
sigmatheta=(-(kappa_y*cos(theta)) + kappa_x*sin(theta))*(K(6,1)/r + power(r,-1 + mVec(1,1))*K(1,1)*(1 + mVec(1,1)) + power(r,-1 + mVec(2,1))*K(2,1)*(1 + mVec(2,1)) + power(r,-1 + mVec(3,1))*K(3,1)*(1 + mVec(3,1)) + power(r,-1 + mVec(4,1))*K(4,1)*(1 + mVec(4,1)) + 3*r*muVec(1,1));
      

      
taurtheta=-(((kappa_x*cos(theta) + kappa_y*sin(theta))*(K(5,1) + K(6,1) + K(6,1)*log(r) + (power(r,mVec(1,1))*K(1,1)*(1 + mVec(1,1)))/mVec(1,1) + (power(r,mVec(2,1))*K(2,1)*(1 + mVec(2,1)))/mVec(2,1) + (power(r,mVec(3,1))*K(3,1)*(1 + mVec(3,1)))/mVec(3,1) + ....
          (power(r,mVec(4,1))*K(4,1)*(1 + mVec(4,1)))/mVec(4,1) + (3*power(r,2)*muVec(1,1))/2.))/r) + ((kappa_x*cos(theta) + kappa_y*sin(theta))*....
      (r*K(5,1) + r*K(6,1)*log(r) + (power(r,1 + mVec(1,1))*K(1,1))/mVec(1,1) + (power(r,1 + mVec(2,1))*K(2,1))/mVec(2,1) + (power(r,1 + mVec(3,1))*K(3,1))/mVec(3,1) + (power(r,1 + mVec(4,1))*K(4,1))/mVec(4,1) + (power(r,3)*muVec(1,1))/2.))/power(r,2);      
      
      
taurz=((kappa_x*cos(theta) + kappa_y*sin(theta))*(power(r,mVec(1,1))*gVec(1,1)*K(1,1) + power(r,mVec(2,1))*gVec(2,1)*K(2,1) + power(r,mVec(3,1))*gVec(3,1)*K(3,1) + power(r,mVec(4,1))*gVec(4,1)*K(4,1) + (K(6,1)*beta(5,6))/beta(6,5) + power(r,2)*muVec(2,1)))/r;


tauthetaz=-((-(kappa_y*cos(theta)) + kappa_x*sin(theta))*(power(r,-1 + mVec(1,1))*gVec(1,1)*K(1,1)*mVec(1,1) + power(r,-1 + mVec(2,1))*gVec(2,1)*K(2,1)*mVec(2,1) + power(r,-1 + mVec(3,1))*gVec(3,1)*K(3,1)*mVec(3,1) + power(r,-1 + mVec(4,1))*gVec(4,1)*K(4,1)*mVec(4,1) + 2*r*muVec(2,1)));



end
