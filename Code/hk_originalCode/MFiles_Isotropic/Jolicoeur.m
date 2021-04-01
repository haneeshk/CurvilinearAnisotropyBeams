% Compliance matrix for isotropic material
clc
close all
clear all



%Geometry
radiiVec=[2e-3; 8e-3; 14e-3;];    % Sequence of Radii
N=length(radiiVec)-1;  % Number of Cylinders

%Imposed displacements
cV=0.0;
epsilon=0.0;

%Material Property
E=7e10;
nu=0.17;

% Set up the material property materices
% CMat=ComposeCMat(N);
CMat=ComposeCMat_Isotropic(N,E,nu);


[betaMat, mMat, muMat, gMat,mpMat, gpMat, UpMat, UppMat, VpMat, WpMat]=SetUpMatPropMat(CMat,N);

% [A, B]= SetUpEquations(mMat, muMat, gMat,mpMat, gpMat, UpMat, UppMat, VpMat, WpMat,....
%                        N, radiiVec, cV, epsilon);


[A, B]= SetUpEquations_NoSlip(mMat, muMat, gMat,mpMat, gpMat, UpMat, UppMat, VpMat, WpMat,....
                       N, radiiVec, cV, epsilon);



X=A\B;
KMat=ComposeKMat(X,N);



EI=ComputeBendingStiffness(KMat, CMat,mMat,gMat,radiiVec,muMat,N);       





rVec=linspace(radiiVec(1,1)+1e-4,radiiVec(N+1,1),1000)';
%rVec=linspace(radiiVec(1,1)+1e-4,radiiVec(N+1,1),1000)';


Mx=10;
theta=pi/2;
kappa_x=Mx/EI;
kappa_y=0;

sigmaz=zeros(length(rVec),1);

sigmar= zeros(length(rVec),1);
sigmatheta=zeros(length(rVec),1);
taurtheta=zeros(length(rVec),1);
taurz=zeros(length(rVec),1);
tauthetaz=zeros(length(rVec),1);


Kn=zeros(6,1);
for i=1:length(rVec)
   r=rVec(i,1);
   
   
       n=find(floor(radiiVec-r)==-1,1,'last'); 
       
       
       
       betan=betaMat(:,:,n);
       Cn=CMat(:,:,n);
       mVecn=mMat(:,n); 
       muVecn=muMat(:,n);
       gVecn=gMat(:,n);
       
       Kn=KMat(:,n);
       Kn(5,1)=0;
       Kn(6,1)=0;

       
       

[sigmaz(i,1), sigmar(i,1), sigmatheta(i,1), taurtheta(i,1), taurz(i,1), tauthetaz(i,1)]=ComputeStress(r,theta,kappa_x,kappa_y,epsilon,betan,Cn,mVecn, muVecn, gVecn,Kn);

end

plot(rVec*1e3,sigmaz*1e-6,'--ks','markers',2)
xlabel('radial position (mm)');
ylabel('axial stress \sigmazz (MPa)');
grid on

figure
plot(rVec*1e3,sigmar*1e-6,'-.k+','markers',2)
xlabel('radial position (mm)');
ylabel('radialstress \sigmarr (MPa)');
grid on
hold all
plot(rVec*1e3,sigmar*1e-6,'--rs','markers',2)



figure
plot(rVec*1e3,taurz*1e-6,'--ks','markers',2)
xlabel('radial position (mm)');
ylabel('radialstress \taurz (MPa)');
grid on


