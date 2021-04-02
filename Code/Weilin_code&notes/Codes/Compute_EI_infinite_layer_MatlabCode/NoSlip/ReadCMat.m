function CMat=ReadCMat(N)
    
fid10=fopen('CMat_Orthotropic2.dat','r+');
CMatn=importdata('CMat_Orthotropic2.dat');
CMat=zeros(6,6,N);

for i=1:N
    CMat(1,1,i)=CMatn.data(1,i);
    CMat(1,2,i)=CMatn.data(2,i);
    CMat(1,3,i)=CMatn.data(3,i);
    CMat(1,4,i)=CMatn.data(4,i);
    CMat(2,2,i)=CMatn.data(5,i);
    CMat(2,3,i)=CMatn.data(6,i);
    CMat(2,4,i)=CMatn.data(7,i);
    CMat(3,3,i)=CMatn.data(8,i);
    CMat(3,4,i)=CMatn.data(9,i);
    CMat(4,4,i)=CMatn.data(10,i);
    CMat(5,5,i)=CMatn.data(11,i);
    CMat(5,6,i)=CMatn.data(12,i);
    CMat(6,6,i)=CMatn.data(13,i);

    CMat(2,1,i)=CMat(1,2,i);
    CMat(3,1,i)=CMat(1,3,i);
    CMat(4,1,i)=CMat(1,4,i);
    CMat(3,2,i)=CMat(2,3,i);
    CMat(4,2,i)=CMat(2,4,i);
    CMat(4,3,i)=CMat(3,4,i);
    CMat(6,5,i)=CMat(5,6,i);
    
    CMat(:,:,i)=CMat(:,:,i)*1e-6;
end

fclose(fid10);

end
