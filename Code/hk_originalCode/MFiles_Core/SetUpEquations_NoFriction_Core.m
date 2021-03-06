function [A, B]= SetUpEquations_NoFriction_Core(mMat, muMat, gMat,mpMat, gpMat, UpMat, UppMat, VpMat, WpMat,N, radiiVec, cV, epsilon);


%Set up system of equations
% AX=B;
A=zeros(7*N, 7*N);
B=zeros(7*N,1);


I=0;  % Equatsion number
%Continuity of \sigma_r
% Set up equation 21
for n=1:N-1
   I=I+1;
      b=radiiVec(1+n,1); 
      B(I,1)=(muMat(3,n+1)-muMat(3,n))*cV*b+(muMat(5,n+1)-muMat(5,n))*epsilon;
    for i=1:2
    J=4*N+2*(n-1)+i;
    A(I,J)=b^(mpMat(i,n)-1);
    
    J=4*N+2*(n+1-1)+i;
    A(I,J)=-b^(mpMat(i,n+1)-1);
    end
end
        


%Continuity of ur, eq. 23, 24
% Set up equation 23
 
for n=1:N-1
   I=I+1;
   b=radiiVec(1+n,1); 
   B(I,1)=(UpMat(5,n+1)-UpMat(5,n))*b^2;
   
   
   J=6*N+n;
   A(I,J)=1;
   
   J=6*N+n+1;
   A(I,J)=-1;
   
   
   for i=1:4
   J=4*(n-1)+i;
   A(I,J)=UpMat(i,n)*b^(mMat(i,n));
   J=4*(n)+i;
   A(I,J)=-UpMat(i,n+1)*b^(mMat(i,n+1));
   
    end
end

% Set up equation 24
for n=1:N-1
   I=I+1;
   b=radiiVec(1+n,1); 
   B(I,1)=(UppMat(3,n+1)-UppMat(3,n))*cV*b^2+(UppMat(4,n+1)-UppMat(4,n))*epsilon*b;
   
   
   
   
   for i=1:2
   J=4*N+2*(n-1)+i;
   A(I,J)=UppMat(i,n)*b^(mpMat(i,n));
      
   J=4*N+2*(n+1-1)+i;
   A(I,J)=-UppMat(i,n+1)*b^(mpMat(i,n+1));
    end
end



% \tau_{r\theta}=0, eqs. 27, 28
% Set up equation 27
 for n=1:N
   I=I+1;
   b=radiiVec(1+n,1); 
   B(I,1)=-muMat(1,n)*b;  

   for i=1:4
  
   J=4*(n-1)+i;
   A(I,J)=b^(mMat(i,n)-1);     
    end
 end
 
 
 
 
% Set up equation 28
for n=1:N-1
I=I+1;
b=radiiVec(1+n,1); 
   B(I,1)=-muMat(1,n+1)*b;  
    
   for i=1:4 
     

       J=4*(n+1-1)+i;
       A(I,J)=b^(mMat(i,n+1)-1);     
    
   end
 end
 
 
 
 % \tau_{rz}=0, eqs. 29, 30
 
% Set up equation 29
for n=1:N
I=I+1;
b=radiiVec(1+n,1); 
B(I,1)=-muMat(2,n)*b;  
    
       for i=1:4
            
            J=4*(n-1)+i;
            A(I,J)=gMat(i,n)*b^(mMat(i,n)-1);     
       end
end


% Set up equation 30
for n=1:N-1
I=I+1;
b=radiiVec(1+n,1); 
B(I,1)=-muMat(2,n+1)*b;  
    for i=1:4
         
         J=4*(n+1-1)+i;
        A(I,J)=gMat(i,n+1)*b^(mMat(i,n+1)-1);     
    
    end
end


% Set up equation 31
n=N;
b=radiiVec(n+1,1);
I=I+1;
B(I,1)=-muMat(3,n)*cV*b-muMat(5,n)*epsilon;
for i=1:2

J=4*N+2*(n-1)+i;
A(I,J)=b^(mpMat(i,n)-1);
end





% cv1=0;
I=I+1;
B(I,1)=0;
J=6*N+1;
A(I,J)=1;


%K(3,1)=0;
I=I+1;
B(I,1)=0;
J=3;
A(I,J)=1;

%K(4,1)=0;
I=I+1;
B(I,1)=0;
J=4;
A(I,J)=1;

%Kp(2,1)=0;
I=I+1;
B(I,1)=0;
J=4*N+2;
A(I,J)=1;



end