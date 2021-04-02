function [A, B]= SetUpEquationsNoSlip(mMat, muMat, gMat,mpMat, gpMat, UpMat, UppMat, VpMat, WpMat,N, radiiVec, cV, epsilon)

% Set up system of equations
% AX=B;
A=zeros(7*N,7*N);
B=zeros(7*N,1);    % (K(1,1),...,K(4,1),...,K(1,N),...,K(4,N),Kp(1,1),Kp(2,1),...,Kp(1,N),Kp(2,N),nu(1),...,nu(N))

I=0;  % Equation number

% Continuity of sigma_r, eq. 20, 21
% Set up equation 20
for n=1:N-1
    I=I+1;
    b=radiiVec(1+n,1);
    B(I,1)=(muMat(1,n+1)-muMat(1,n))*b;
    for i=1:4
        J=4*(n-1)+i;
        A(I,J)=b^(mMat(i,n)-1);
        J=4*(n-1+1)+i;
        A(I,J)=-b^(mMat(i,n+1)-1);
    end
end

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


% Continuity of tau_rtheta
% Set up equation 22
for n=1:N-1
    I=I+1;
    b=radiiVec(1+n,1);
    B(I,1)=(muMat(2,n+1)-muMat(2,n))*b;
    for i=1:4
        J=4*(n-1)+i;
        A(I,J)=gMat(i,n)*b^(mMat(i,n)-1);
        J=4*(n-1+1)+i;
        A(I,J)=-gMat(i,n+1)*b^(mMat(i,n+1)-1);
    end
end

% Continuity of u_r, eq. 23, 24
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
        J=4*(n-1+1)+i;
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


% Continuity of u_theta, eq. 25
% Set up equation 25
for n=1:N-1
    I=I+1;
    b=radiiVec(1+n,1); 
    B(I,1)=(VpMat(5,n+1)-VpMat(5,n))*b^2;

    J=6*N+n;
    A(I,J)=1;
    J=6*N+n+1;
    A(I,J)=-1;
    
    for i=1:4
        J=4*(n-1)+i;
        A(I,J)=VpMat(i,n)*b^(mMat(i,n));
        J=4*(n-1+1)+i;
        A(I,J)=-VpMat(i,n+1)*b^(mMat(i,n+1));
    end
end

% Continuity of u_theta, eq. 26
% Set up equation 26
for n=1:N-1
    I=I+1;
    b=radiiVec(1+n,1); 
    B(I,1)=(WpMat(5,n+1)-WpMat(5,n))*b^2;
    
    for i=1:4
        J=4*(n-1)+i;
        A(I,J)=WpMat(i,n)*b^(mMat(i,n));
        J=4*(n-1+1)+i;
        A(I,J)=-WpMat(i,n+1)*b^(mMat(i,n+1));
    end
    
end


% free inner and outer surfaces, n=0 and N
% n=N;
n=N;
b=radiiVec(n+1,1);

% Set up equation 27
I=I+1;
B(I,1)=-muMat(1,n)*b;  
for i=1:4
    J=4*(n-1)+i;
    A(I,J)=b^(mMat(i,n)-1);
end

% Set up equation 29
I=I+1;
B(I,1)=-muMat(2,n)*b;  
for i=1:4
    J=4*(n-1)+i;
    A(I,J)=gMat(i,n)*b^(mMat(i,n)-1);     
end

% Set up equation 31
I=I+1;
B(I,1)=-muMat(3,n)*cV*b-muMat(5,n)*epsilon;
for i=1:2
    J=4*N+2*(n-1)+i;
    A(I,J)=b^(mpMat(i,n)-1);
end

%n=0
n=0;
b=radiiVec(n+1,1);

% Set up equation 28
I=I+1;
B(I,1)=-muMat(1,n+1)*b;  
for i=1:4
    J=4*(n-1+1)+i;
    A(I,J)=b^(mMat(i,n+1)-1);
end

% Set up equation 30
I=I+1;
B(I,1)=-muMat(2,n+1)*b;  
for i=1:4
    J=4*(n-1+1)+i;
    A(I,J)=gMat(i,n+1)*b^(mMat(i,n+1)-1);
end

% Set up equation 32
I=I+1;
B(I,1)=-muMat(3,n+1)*cV*b-muMat(5,n+1)*epsilon;
for i=1:2
    J=4*N+2*(n+1-1)+i;
    A(I,J)=b^(mpMat(i,n+1)-1);
end

% cV1=0;
I=I+1;
B(I,1)=0;
J=6*N+1;
A(I,J)=1;

end