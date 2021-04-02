% Compliance matrix for isotropic material
%clc
close all
clear all



%Geometry
%radiiVec=[0;  8e-3; 14e-3];    % Sequence of Radii
radiiVec=[0; 14e-3;]; 

%radiiVec=linspace(0,14e-3,7);
%radiiVec=linspace(0,14e-3,9);
%radiiVec=linspace(0,14e-3,11);
%  radiiVec=linspace(0,14e-3,8);
%  radiiVec=radiiVec';


N=length(radiiVec)-1;  % Number of Cylinders

%Imposed displacements
cV=0.0;
epsilon=0.0;




% Set up the material property materices
%CMat=ComposeCMat(N); 
CMat=ComposeCMat_Homo(N); 


[betaMat, mMat, muMat, gMat,mpMat, gpMat, UpMat, UppMat, VpMat, WpMat]=SetUpMatPropMat(CMat,N);

%!!!Note that m and mp value may cause probelm for the core at r=0. In this
%example, m(3,1)+2 is accidentally very close to 0, and this caused a
%problem of 0*Na
%while calculation the EI. Hence, the value of EI is set to 0 manually for
%the calculation at that two points. m+1 does cause similar
%problem during later calculation of stress, but it only fails the first
%point at r=0, without affecting later values.






  [A, B]= SetUpEquations_NoFriction_Core(mMat, muMat, gMat,mpMat, gpMat, UpMat, UppMat, VpMat, WpMat,....
                         N, radiiVec, cV, epsilon);


%    [A, B]= SetUpEquations_NoSlip_Core(mMat, muMat, gMat,mpMat, gpMat, UpMat, UppMat, VpMat, WpMat,....
%                           N, radiiVec, cV, epsilon);

                        
%  [A, B]= SetUpEquations(mMat, muMat, gMat,mpMat, gpMat, UpMat, UppMat, VpMat, WpMat,....
%                         N, radiiVec, cV, epsilon);



X=A\B;
%X=pcg(A,B);
%X=bicg(A,B);
%X=bicgstab(A,B);
%X=cgs(A,B);
%X=gmres(A,B,10);
%X=lsqr(A,B);
%X=minres(A,B);
%X=qmr(A,B);
%X=symmlq(A,B);
%X=linsolve(A,B);


KMat=ComposeKMat(X,N);
KpMat=ComposeKpMat(X,N);
nuMat=ComposenuMat(X,N);


EI=ComputeBendingStiffness(KMat, CMat,mMat,gMat,radiiVec,muMat,N);       





rVec=linspace(radiiVec(1,1),radiiVec(N+1,1),10000)';
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
   
        
        for m=1:N
            
            if r<=radiiVec(m+1,1)
            n=m;
            
            break;
            end
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
