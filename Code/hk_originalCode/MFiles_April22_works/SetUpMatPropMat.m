function [ betaMat, mMat, muMat, gMat,mpMat, gpMat, UpMat, UppMat, VpMat, WpMat]=SetUpMatPropMat(CMat,N)

% Initialization
 
 betaMat=zeros(6,6,N);
 mMat=zeros(4,N);
 muMat=zeros(5,N);
 gMat=zeros(4,N);
 mpMat=zeros(2,N);
 gpMat=zeros(2,N);
 UpMat=zeros(5,N);
 UppMat=zeros(4,N);
 VpMat=zeros(5,N);
 WpMat=zeros(5,N);
 




for n=1:N


Cn=CMat(:,:,n) % Compliance Tensor for Cylinder n
    

[ betan, mVecn, muVecn, gVecn, mpVecn, gpVecn, UpVecn, UppVecn, VpVecn, WpVecn]= ComputeMaterialProperties(Cn);

betaMat(:,:,n)=betan;
mMat(:,n)=mVecn; 
muMat(:,n)=muVecn; 
gMat(:,n)=gVecn;
mpMat(:,n)=mpVecn;
gpMat(:,n)=gpVecn;
UpMat(:,n)=UpVecn;
UppMat(:,n)=UppVecn;
VpMat(:,n)=VpVecn;
WpMat(:,n)=WpVecn;
clear Cn, betan, mVecn, muVecn, gVecn, mpVecn, gpVecn, UpVecn, UppVecn, VpVecn, WpVecn;
end

end







function [ beta, mVec, muVec, gVec, mpVec, gpVec, UpVec, UppVec, VpVec, WpVec]= ComputeMaterialProperties(C)
nDim=6;




% equation 2 
beta=zeros(nDim,nDim)
for i=1:nDim
    for j=1:nDim
        beta(i,j)=C(i,j)-(C(i,3)*C(3,j)/C(3,3))
    end
end





% equation 12
mVec=zeros(4,1);


a=beta(2,2)*beta(4,4)-beta(2,4)^2;
b=beta(2,4)*(2*beta(1,4)+beta(2,4)+2*beta(5,6))....
   -beta(4,4)*(beta(1,1)+2*beta(1,2)+beta(2,2)+beta(6,6))...
   -beta(2,2)*beta(5,5)....
   +beta(1,4)^2;
c=beta(5,5)*(beta(1,1)+2*beta(1,2)+beta(2,2)+beta(6,6))-beta(5,6)^2




mVec(1,1)=sqrt(....
            (-b+sqrt(b^2-4*a*c))/(2*a)....
        );
    
    
mVec(2,1)=sqrt(....
            (-b-sqrt(b^2-4*a*c))/(2*a)....
        );

mVec(3,1)=-sqrt(....
            (-b+sqrt(b^2-4*a*c))/(2*a)....
        );
    
    
mVec(4,1)=-sqrt(....
            (-b-sqrt(b^2-4*a*c))/(2*a)....
        );


    
muVec=zeros(5,1);
    
% equation 13,
BetaMat12=zeros(2,2);
BetaMat12(1,1)=-2*beta(1,4)-6*beta(2,4)+beta(5,6);
BetaMat12(1,2)=4*beta(4,4)-beta(5,5);
BetaMat12(2,1)=-beta(1,1)-2*beta(1,2)+3*beta(2,2)-beta(6,6);
BetaMat12(2,2)=2*beta(1,4)-2*beta(2,4)+beta(5,6);

mu12=inv(BetaMat12)*(1/C(3,3))*[2*C(3,4); C(1,3)-C(2,3)]



BetaMat34=zeros(2,2);
BetaMat34(1,1)=beta(1,4)+2*beta(2,4);
BetaMat34(1,2)=-beta(4,4);
BetaMat34(2,1)=4*beta(2,2)-beta(1,1);
BetaMat34(2,2)=beta(1,4)-2*beta(2,4);


mu34=inv(BetaMat34)*[1; 0];
%mu5=0;
mu5=C(3,4)*(beta(2,4)-beta(1,4))+beta(4,4)*(C(1,3)-C(2,3))
 mu5=mu5/(....
 C(3,3)*....
 (....
 beta(1,4)^2-beta(2,4)^2+beta(4,4)*(beta(2,2)-beta(1,1))....
 )....
 );



muVec=[mu12; mu34; mu5]


% Computing g
gVec=zeros(4,1)
 for i=1:4
  gVec(i,1)=beta(2,4)*mVec(i,1)^2+(beta(1,4)+beta(2,4))*mVec(i,1)-beta(5,6);
  gVec(i,1)=gVec(i,1)/(beta(4,4)*mVec(i,1)^2-beta(5,5));
 end


mpVec=zeros(2,1);
mpVec(1,1)=(beta(1,1)*beta(4,4)-beta(1,4)^2);
mpVec(1,1)=mpVec(1,1)/(beta(2,2)*beta(4,4)-beta(2,4)^2);
mpVec(1,1)=abs(sqrt(mpVec(1,1)));
mpVec(2,1)=-mpVec(1,1);



gpVec=zeros(2,1);
 for i=1:2
  gpVec(i,1)=(beta(1,4)+beta(2,4)*mpVec(i,1))/beta(4,4);
 end




 UpVec=zeros(5,1);
 for i=1:4
 UpVec(i,1)=beta(1,1)+beta(1,2)*(mVec(i,1)+1)-beta(1,4)*gVec(i,1)*mVec(i,1);
 UpVec(i,1)=UpVec(i,1)/mVec(i,1);
 end
 
 UpVec(5,1)=muVec(1,1)*(beta(1,1)+3*beta(1,2));
 UpVec(5,1)=UpVec(5,1)-2*beta(1,4)*muVec(2,1);
 UpVec(5,1)=UpVec(5,1)+C(1,3)/C(3,3);
 UpVec(5,1)= UpVec(5,1)/2;
 
 
 
 UppVec=zeros(4,1);
  for i=1:2
UppVec(i,1)=beta(1,1);
 UppVec(i,1)=UppVec(i,1)+beta(1,2)*mpVec(i,1);
 UppVec(i,1)=UppVec(i,1)-beta(1,4)*gpVec(i,1);
 UppVec(i,1)=UppVec(i,1)/mpVec(i,1);
 end
 
 UppVec(3,1)=muVec(3,1)*(beta(1,1)+2*beta(1,2));
 UppVec(3,1)=UppVec(3,1)-muVec(4,1)*beta(1,4);
 UppVec(3,1)=UppVec(3,1)/2;
 
 
 
UppVec(4,1)=(C(1,3)-(beta(1,4)/beta(4,4))*C(3,4));
UppVec(4,1)=UppVec(4,1)/C(3,3)
UppVec(4,1)=UppVec(4,1)+muVec(5,1)*(beta(1,1)+beta(1,2)-beta(1,4)*((beta(1,4)+beta(2,4))/(beta(4,4))));
 




VpVec=zeros(5,1)
for i=1:4
VpVec(i,1)=(beta(1,1)+beta(1,2)-beta(2,2)*mVec(i,1)*(mVec(i,1)+1)-gVec(i,1)*mVec(i,1)*(beta(1,4)-beta(2,4)*mVec(i,1)));
VpVec(i,1)=VpVec(i,1)/mVec(i,1);
end

VpVec(5,1)=(muVec(1,1)*(beta(1,1)+beta(1,2)-6*beta(2,2))-2*muVec(2,1)*(beta(1,4)-2*beta(2,4))+((C(1,3)-2*C(2,3))/C(3,3)))/2;




WpVec=zeros(5,1);
for i=1:4
WpVec(i,1)=(beta(5,5)*gVec(i,1)-beta(5,6))/mVec(i,1);
end

WpVec(5,1)=(beta(5,5)*muVec(2,1)-beta(5,6)*muVec(1,1))/2;

end