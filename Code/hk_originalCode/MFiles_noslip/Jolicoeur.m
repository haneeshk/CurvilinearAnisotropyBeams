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




% Set up the material property materices
CMat=ComposeCMat(N);

[betaMat, mMat, muMat, gMat,mpMat, gpMat, UpMat, UppMat, VpMat, WpMat]=SetUpMatPropMat(CMat,N);

% [A, B]= SetUpEquations(mMat, muMat, gMat,mpMat, gpMat, UpMat, UppMat, VpMat, WpMat,....
%                        N, radiiVec, cV, epsilon);


[A, B]= SetUpEquations_NoSlip(mMat, muMat, gMat,mpMat, gpMat, UpMat, UppMat, VpMat, WpMat,....
                       N, radiiVec, cV, epsilon);



X=A\B;
KMat=ComposeKMat(X,N);
KpMat=ComposeKpMat(X,N);
nuMat=ComposenuMat(X,N);


EI=ComputeBendingStiffness(KMat, CMat,mMat,gMat,radiiVec,muMat,N);       





rVec=linspace(radiiVec(1,1),radiiVec(N+1,1),1000)';
%rVec=linspace(radiiVec(1,1)+1e-4,radiiVec(N+1,1),1000)';


Mx=10;
theta=0*pi/2;
kappa_x=Mx/EI;
kappa_y=0;

sigmaz=zeros(length(rVec),1);

sigmar= zeros(length(rVec),1);
sigmatheta=zeros(length(rVec),1);
taurtheta=zeros(length(rVec),1);
taurz=zeros(length(rVec),1);
tauthetaz=zeros(length(rVec),1);
U=zeros(length(rVec),1);
V=zeros(length(rVec),1);
W=zeros(length(rVec),1);
ur=zeros(length(rVec),1);
utheta=zeros(length(rVec),1);
w=zeros(length(rVec),1);


Kn=zeros(6,1);
Kpn=zeros(2,1);
nun=zeros(1,1);

z=0; %position for displacement

for i=1:length(rVec)
   r=rVec(i,1);
   
        if r==radiiVec(1,1)
            n=1;
        else
            n=find(floor(radiiVec-r)==-1,1,'last');
        end
       
       
       betan=betaMat(:,:,n);
       Cn=CMat(:,:,n);
       mVecn=mMat(:,n); 
       muVecn=muMat(:,n);
       gVecn=gMat(:,n);
       
       Kn=KMat(:,n);
       Kn(5,1)=0;
       Kn(6,1)=0;
       
       Kpn=KpMat(:,n);
       
       nun=nuMat(n,1);

       Upn=UpMat(:,n);
       Uppn=UppMat(:,n);
       Vpn=VpMat(:,n);
       Wpn=WpMat(:,n);
       

[sigmaz(i,1), sigmar(i,1), sigmatheta(i,1), taurtheta(i,1), taurz(i,1), tauthetaz(i,1)]=ComputeStress(r,theta,kappa_x,kappa_y,epsilon,betan,Cn,mVecn, muVecn, gVecn,Kn);


[U(i,1),V(i,1),W(i,1),ur(i,1),utheta(i,1),w(i,1)]=ComputeDisplacement(kappa_x,kappa_y,theta,r,Kn,Kpn,cV,epsilon,Upn,Uppn,Vpn,Wpn,mVecn,muVecn,z,nun);



end

PlotResult(rVec,sigmaz,sigmar,taurz,taurtheta,ur,utheta,w)

% figure(4)
% hold all
% plot(rVec*1e3,sigmar*1e-6,'-.k+','markers',2)
% xlabel('Radial Position (mm)');
% ylabel('Radial Stress \sigma_{rr} (MPa)');
% grid on
