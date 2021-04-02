function CMat=ComposeCMat_Isotropic(N,E,nu)



nDim=6;
C=zeros(nDim,nDim);  % Compliance tensor
C(1,1)=1;
C(1,2)=-nu;
C(1,3)=-nu;

C(2,1)=-nu;
C(2,2)=1;
C(2,3)=-nu;


C(3,1)=-nu;
C(3,2)=-nu;
C(3,3)=1;



C(4,4)=2*(1+nu);
C(5,5)=2*(1+nu);
C(6,6)=2*(1+nu);



C=C/E;

CMat=zeros(6,6,N);

for n=1:N

    CMat(:,:,n)=C;

end


end