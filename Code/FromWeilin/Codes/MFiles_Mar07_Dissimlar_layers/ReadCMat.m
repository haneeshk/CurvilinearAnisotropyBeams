function CMat=ReadCMat(N,alpha,flag_mat,flag_sym)

CMat=zeros(6,6,N);
CMat0=zeros(6,6,1);
% Rotation matrix about Er axis
R=zeros(6,6); s=sin(alpha);c=cos(alpha);
R(1,1)=1; R(1,2:6)=0;
R(2,1)=0; R(2,5:6)=0; R(2,2)=c^2; R(2,3)=s^2; R(2,4)=2*c*s;
R(3,1)=0; R(3,5:6)=0; R(3,2)=s^2; R(3,3)=c^2; R(3,4)=-2*c*s;
R(4,1)=0; R(4,5:6)=0; R(4,2)=-c*s; R(4,3)=c*s; R(4,4)=c^2-s^2;
R(5,1:4)=0; R(5,5)=c; R(5,6)=-s;
R(6,1:4)=0; R(6,5)=s; R(6,6)=c;

Rp=zeros(6,6); sp=sin(alpha+40/180*pi);cp=cos(alpha+40/180*pi);
Rp(1,1)=1; Rp(1,2:6)=0;
Rp(2,1)=0; Rp(2,5:6)=0; Rp(2,2)=cp^2; Rp(2,3)=sp^2; Rp(2,4)=2*cp*sp;
Rp(3,1)=0; Rp(3,5:6)=0; Rp(3,2)=sp^2; Rp(3,3)=cp^2; Rp(3,4)=-2*cp*sp;
Rp(4,1)=0; Rp(4,5:6)=0; Rp(4,2)=-cp*sp; Rp(4,3)=cp*sp; Rp(4,4)=cp^2-sp^2;
Rp(5,1:4)=0; Rp(5,5)=cp; Rp(5,6)=-sp;
Rp(6,1:4)=0; Rp(6,5)=sp; Rp(6,6)=cp;

% fid10=fopen('CMat_CurvOrtho105.dat','r+');
% CMatn=importdata('CMat_CurvOrtho105.dat');

% fid10=fopen('CMat_CurvOrtho90.dat','r+');
% CMatn=importdata('CMat_CurvOrtho90.dat');

% fid10=fopen('CMat_TransIso.dat','r+');
% CMatn=importdata('CMat_TransIso.dat');

% fid10=fopen('CMat_Isotropic.dat','r+');
% CMatn=importdata('CMat_Isotropic.dat');

if flag_mat==2 % heterogeneous layers
%     fid10000=fopen('CMat_Orthotropic2.dat','r+');
%     CMatn=importdata('CMat_Orthotropic2.dat');

    fid10000=fopen('CMat_CurvOrtho.dat','r+');
    CMatn=importdata('CMat_CurvOrtho.dat');
for i=1:1
    ii=1;
    CMat0(1,1,i)=CMatn.data(1,ii);
    CMat0(1,2,i)=CMatn.data(2,ii);
    CMat0(1,3,i)=CMatn.data(3,ii);
    CMat0(1,4,i)=CMatn.data(4,ii);
    CMat0(2,2,i)=CMatn.data(5,ii);
    CMat0(2,3,i)=CMatn.data(6,ii);
    CMat0(2,4,i)=CMatn.data(7,ii);
    CMat0(3,3,i)=CMatn.data(8,ii);
    CMat0(3,4,i)=CMatn.data(9,ii);
    CMat0(4,4,i)=CMatn.data(10,ii);
    CMat0(5,5,i)=CMatn.data(11,ii);
    CMat0(5,6,i)=CMatn.data(12,ii);
    CMat0(6,6,i)=CMatn.data(13,ii);

    CMat0(2,1,i)=CMat0(1,2,i);
    CMat0(3,1,i)=CMat0(1,3,i);
    CMat0(4,1,i)=CMat0(1,4,i);
    CMat0(3,2,i)=CMat0(2,3,i);
    CMat0(4,2,i)=CMat0(2,4,i);
    CMat0(4,3,i)=CMat0(3,4,i);
    CMat0(6,5,i)=CMat0(5,6,i);

%     CMat0(:,:,i)=CMat0(:,:,i)*1e-6;
end
    
    for i=1:N
        ii=mod(i,2); % two different materials: ii=1 - mat1; ii=2 - mat2
        if ii==1
            % stiffness matrix in cylindrical coordinates {Er,Etheta,Ez}
            Ri=inv(R); Rit=transpose(Ri);
            CMat(:,:,i)=Rit*CMat0(:,:,1)*Ri;
        elseif ii==0
            Rpi=inv(Rp); Rpit=transpose(Rpi);
            CMat(:,:,i)=Rpit*CMat0(:,:,1)*Rpi;
        end        
    end
    
    fclose(fid10000);
    
elseif flag_mat==1 % homogeneous layers
    if flag_sym==1 % orthotropic
        fid20000=fopen('CMat_CurvOrtho.dat','r+');
        CMatn=importdata('CMat_CurvOrtho.dat');
        fclose(fid20000);
    elseif flag_sym==2 % transverse isotropic
        fid30000=fopen('CMat_TransIso.dat','r+');
        CMatn=importdata('CMat_TransIso.dat');
        fclose(fid30000);
    elseif flag_sym==3 % isotropic
        fid40000=fopen('CMat_Isotropic.dat','r+');
        CMatn=importdata('CMat_Isotropic.dat');
        fclose(fid40000);
    end
    for i=1:1
        ii=1;
        CMat(1,1,i)=CMatn.data(1,ii);
        CMat(1,2,i)=CMatn.data(2,ii);
        CMat(1,3,i)=CMatn.data(3,ii);
        CMat(1,4,i)=CMatn.data(4,ii);
        CMat(2,2,i)=CMatn.data(5,ii);
        CMat(2,3,i)=CMatn.data(6,ii);
        CMat(2,4,i)=CMatn.data(7,ii);
        CMat(3,3,i)=CMatn.data(8,ii);
        CMat(3,4,i)=CMatn.data(9,ii);
        CMat(4,4,i)=CMatn.data(10,ii);
        CMat(5,5,i)=CMatn.data(11,ii);
        CMat(5,6,i)=CMatn.data(12,ii);
        CMat(6,6,i)=CMatn.data(13,ii);

        CMat(2,1,i)=CMat(1,2,i);
        CMat(3,1,i)=CMat(1,3,i);
        CMat(4,1,i)=CMat(1,4,i);
        CMat(3,2,i)=CMat(2,3,i);
        CMat(4,2,i)=CMat(2,4,i);
        CMat(4,3,i)=CMat(3,4,i);
        CMat(6,5,i)=CMat(5,6,i);

%         CMat(:,:,i)=CMat(:,:,i)*1e-6;
    end
    
    % stiffness matrix in cylindrical coordinates {Er,Etheta,Ez}
    Ri=inv(R); Rit=transpose(Ri);
    CMat(:,:,1)=Rit*CMat(:,:,1)*Ri;
    
    for i=2:N
        CMat(:,:,i)=CMat(:,:,1);
    end
    
else
    display('Ilegal input!');
end

end
