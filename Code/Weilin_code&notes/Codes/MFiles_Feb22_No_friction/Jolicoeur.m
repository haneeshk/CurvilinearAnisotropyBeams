% Compliance matrix for isotropic material
% clc
% close all
% clear all
% 
% %Geometry
% radiiVec=[2e-3; 8e-3; 14e-3];    % Sequence of Radii
% N=length(radiiVec)-1;             % Number of Cylinders

function EI = Jolicoeur(radiiVec,I)

N=length(radiiVec)-1;

% Set up the material property materices
% CMat=ComposeCMat(N);
CMat=ReadCMat(N);

%Imposed displacements
cV=0.0;
epsilon=0.0;

[betaMat, mMat, muMat, gMat,mpMat, gpMat, UpMat, UppMat, VpMat, WpMat]=SetUpMatPropMat(CMat,N);

%[A, B]= SetUpEquationsNoSlip(mMat, muMat, gMat,mpMat, gpMat, UpMat, UppMat, VpMat, WpMat,N, radiiVec, cV, epsilon);
[A, B]= SetUpEquationsNoFriction(mMat, muMat, gMat,mpMat, gpMat, UpMat, UppMat, VpMat, WpMat,N, radiiVec, cV, epsilon);


X=A\B;
%X=pinv(A)*B;

KMat=ComposeKMat(X,N);
KpMat=ComposeKpMat(X,N);
nuMat=ComposenuMat(X,N);

EI=ComputeBendingStiffness(KMat,CMat,mMat,gMat,radiiVec,muMat,N);       

rVec=linspace(radiiVec(1,1),radiiVec(N+1,1),1000)';
theta=pi/2;
%theta=0;

Mx=10;
My=0;
kappa_x=Mx/EI;
kappa_y=My/EI;

sigmaz=zeros(length(rVec),1);
sigmar= zeros(length(rVec),1);
sigmatheta=zeros(length(rVec),1);
taurtheta=zeros(length(rVec),1);
taurz=zeros(length(rVec),1);
tauthetaz=zeros(length(rVec),1);

ur=zeros(length(rVec),1);
utheta=zeros(length(rVec),1);
w=zeros(length(rVec),1);

% fid1=fopen('Stress.dat','w');
% fid2=fopen('Displacement.dat','w');
% fid1=fopen('Stress.dat','a+');
% fid2=fopen('Displacement.dat','a+');
% fprintf(fid1,'r         sigmar      sigmatheta   sigmaz  taurtheta    taurz    tauthetaz\r\n');
% fprintf(fid2,'r         ur            utheta            w\r\n');

Kn=zeros(4,1);
Kpn=zeros(2,1);

z=0; %position for displacement

for i=1:length(rVec)
   r=rVec(i,1);

   if r==radiiVec(1,1)
       n=1;
   else
       n=find(floor(radiiVec-r)==-1,1,'last');
   end       
   
   Cn=CMat(:,:,n);
   betan=betaMat(:,:,n);
   mVecn=mMat(:,n); 
   mpVecn=mpMat(:,n);
   muVecn=muMat(:,n);
   gVecn=gMat(:,n);
   gpVecn=gpMat(:,n);
   UpVecn=UpMat(:,n);
   UppVecn=UppMat(:,n);
   VpVecn=VpMat(:,n);
   WpVecn=WpMat(:,n);
   
   Kn=KMat(:,n);
   Kpn=KpMat(:,n);
   nun=nuMat(:,n);
   
   [sigmaz(i,1), sigmar(i,1), sigmatheta(i,1), taurtheta(i,1), taurz(i,1), tauthetaz(i,1)]=ComputeStressNew(r,theta,kappa_x,kappa_y,epsilon,cV,betan,Cn,mVecn,mpVecn,muVecn,gVecn,gpVecn,Kn,Kpn);

   [ur(i,1),utheta(i,1),w(i,1)]=ComputeDisplacement(kappa_x,kappa_y,theta,r,Kn,Kpn,cV,epsilon,UpVecn,UppVecn,VpVecn,WpVecn,mVecn,muVecn,nun,z);
   
%    fprintf(fid1,'%8.5f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f \n',r*1.0E3,sigmar(i,1)*1.0E-6,sigmatheta(i,1)*1.0E-6,sigmaz(i,1)*1.0E-6,...
%                  taurtheta(i,1)*1.0E-6, taurz(i,1)*1.0E-6, tauthetaz(i,1)*1.0E-6);
%    fprintf(fid2,'%8.5f %12.10f %12.10f %12.10f \n',r*1.0E3,ur(i,1)*1.0E3,utheta(i,1)*1.0E3,w(i,1)*1.0E3);         

end

% fclose(fid1);
% fclose(fid2);
% if I==15
%     PlotResult(rVec,sigmar,sigmatheta,sigmaz,taurz,taurtheta,tauthetaz,ur,utheta,w);
% end
end