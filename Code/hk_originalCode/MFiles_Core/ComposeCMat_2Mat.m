function CMat=ComposeCMat_2Mat(N,nu1,E1,nu2,E2)



C1(1,1)=1;
C1(1,2)=-nu1;
C1(1,3)=-nu1;
C1(2,2)=1;
C1(2,3)=-nu1;
C1(3,3)=1;
C1(4,4)=1+nu1;
C1(5,5)=1+nu1;
C1(6,6)=1+nu1;

C1(1,4)=0.00003;
C1(2,4)=0.00002;
C1(3,4)=0.00001;
C1(5,6)=0.00004;



C1=C1/E1;

for i=1:6
    for j=1:6
        C1(j,i)=C1(i,j);
    end
end




C2(1,1)=1;
C2(1,2)=-nu2;
C2(1,3)=-nu2;
C2(2,2)=1;
C2(2,3)=-nu2;
C2(3,3)=1;
C2(4,4)=1+nu2;
C2(5,5)=1+nu2;
C2(6,6)=1+nu2;

C2(1,4)=0.00003;
C2(2,4)=0.00002;
C2(3,4)=0.00001;
C2(5,6)=0.00004;



C2=C2/E2;

for i=1:6
    for j=1:6
        C2(j,i)=C2(i,j);
    end
end


% nu3=0.17;
% E3=20*10^9;
% 
% 
% C3(1,1)=1;
% C3(1,2)=-nu3;
% C3(1,3)=-nu3;
% C3(2,2)=1;
% C3(2,3)=-nu3;
% C3(3,3)=1;
% C3(4,4)=1+nu3;
% C3(5,5)=1+nu3;
% C3(6,6)=1+nu3;
% 
% C3(1,4)=0.00003;
% C3(2,4)=0.00002;
% C3(3,4)=0.00001;
% C3(5,6)=0.00004;
% 
% 
% 
% C3=C3/E3;
% 
% for i=1:6
%     for j=1:6
%         C3(j,i)=C3(i,j);
%     end
% end












% C1(1,1)= 1.052632e-4;
% C1(1,2)=-6.456815e-6;
% C1(1,3)=-8.280027e-6;
% C1(1,4)= 1.052632e-6;
% C1(2,2)= 1.263265e-4;
% C1(2,3)=-2.463235e-5;
% C1(2,4)=-2.831390e-5;
% C1(3,3)= 4.176170e-5;
% C1(3,4)= 7.713743e-5;
% C1(4,4)= 3.391176e-4;
% C1(5,5)= 4.832532e-4;
% C1(5,6)=-6.250000e-5;
% C1(6,6)= 2.667468e-4;
% 
% 
% C1(2,1)=C1(1,2);
% C1(3,1)=C1(1,3);
% C1(4,1)=C1(1,4);
% C1(3,2)=C1(2,3);
% C1(4,2)=C1(2,4);
% C1(4,3)=C1(3,4);
% C1(6,5)=C1(5,6);
% 
% 
% 
% 
% C2(1,1)= 1.052632e-4;
% C2(1,2)=-6.691803e-6;
% C2(1,3)=-8.045040e-6;
% C2(1,4)=-1.612725e-6;
% C2(2,2)= 1.359339e-4;
% C2(2,3)=-4.513900e-5;
% C2(2,4)= 2.255642e-5;
% C2(3,3)= 7.316760e-5;
% C2(3,4)=-9.735841e-5;
% C2(4,4)= 2.570911e-4;
% C2(5,5)= 4.553485e-4;
% C2(5,6)= 9.575556e-5;
% C2(6,6)= 2.946515e-4;
% 
% C2(2,1)=C2(1,2);
% C2(3,1)=C2(1,3);
% C2(4,1)=C2(1,4);
% C2(3,2)=C2(2,3);
% C2(4,2)=C2(2,4);
% C2(4,3)=C2(3,4);
% C2(6,5)=C2(5,6);
% 
% 
% C1=C1*1e-6;
% C2=C2*1e-6;

CMat=zeros(6,6,N);

for n=1:N

    if rem(n,2)==1
        CMat(:,:,n)=C1;
    else
            CMat(:,:,n)=C2;
    end

end


end