function CMat=ReadCMat(N,flag)

% fid10=fopen('CMat_CurvOrtho105.dat','r+');
% CMatn=importdata('CMat_CurvOrtho105.dat');

% fid10=fopen('CMat_CurvOrtho90.dat','r+');
% CMatn=importdata('CMat_CurvOrtho90.dat');

% fid10=fopen('CMat_TransIso.dat','r+');
% CMatn=importdata('CMat_TransIso.dat');

% fid10=fopen('CMat_Isotropic.dat','r+');
% CMatn=importdata('CMat_Isotropic.dat');

fid10=fopen('CMat_Orthotropic2.dat','r+');
CMatn=importdata('CMat_Orthotropic2.dat');

CMat=zeros(6,6,N);
if flag==2
    for i=1:N
        ii=mod(i,2);
        if ii==1
            CMat(1,1,i)=CMatn.data(1,1);
            CMat(1,2,i)=CMatn.data(2,1);
            CMat(1,3,i)=CMatn.data(3,1);
            CMat(1,4,i)=CMatn.data(4,1);
            CMat(2,2,i)=CMatn.data(5,1);
            CMat(2,3,i)=CMatn.data(6,1);
            CMat(2,4,i)=CMatn.data(7,1);
            CMat(3,3,i)=CMatn.data(8,1);
            CMat(3,4,i)=CMatn.data(9,1);
            CMat(4,4,i)=CMatn.data(10,1);
            CMat(5,5,i)=CMatn.data(11,1);
            CMat(5,6,i)=CMatn.data(12,1);
            CMat(6,6,i)=CMatn.data(13,1);
    
            CMat(2,1,i)=CMat(1,2,i);
            CMat(3,1,i)=CMat(1,3,i);
            CMat(4,1,i)=CMat(1,4,i);
            CMat(3,2,i)=CMat(2,3,i);
            CMat(4,2,i)=CMat(2,4,i);
            CMat(4,3,i)=CMat(3,4,i);
            CMat(6,5,i)=CMat(5,6,i);
    
%             CMat(:,:,i)=CMat(:,:,i)*1e-6;
        elseif ii==0
            CMat(1,1,i)=CMatn.data(1,2);
            CMat(1,2,i)=CMatn.data(2,2);
            CMat(1,3,i)=CMatn.data(3,2);
            CMat(1,4,i)=CMatn.data(4,2);
            CMat(2,2,i)=CMatn.data(5,2);
            CMat(2,3,i)=CMatn.data(6,2);
            CMat(2,4,i)=CMatn.data(7,2);
            CMat(3,3,i)=CMatn.data(8,2);
            CMat(3,4,i)=CMatn.data(9,2);
            CMat(4,4,i)=CMatn.data(10,2);
            CMat(5,5,i)=CMatn.data(11,2);
            CMat(5,6,i)=CMatn.data(12,2);
            CMat(6,6,i)=CMatn.data(13,2);
    
            CMat(2,1,i)=CMat(1,2,i);
            CMat(3,1,i)=CMat(1,3,i);
            CMat(4,1,i)=CMat(1,4,i);
            CMat(3,2,i)=CMat(2,3,i);
            CMat(4,2,i)=CMat(2,4,i);
            CMat(4,3,i)=CMat(3,4,i);
            CMat(6,5,i)=CMat(5,6,i);
    
%             CMat(:,:,i)=CMat(:,:,i)*1e-6;
        end
        
    end
    
elseif flag==1
    
    for i=1:1
        ii=1; % same properties, choose in kayer 1 or 2
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

    for i=2:N
        CMat(:,:,i)=CMat(:,:,1);
    end
    
else
    display('Ilegal input!');
    exit;
end

fclose(fid10);

end
