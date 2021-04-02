function nuMat=ComposenuMat(X,N)

nuMat=zeros(N,1);

for n=1:N
        J=6*N+n;
        nuMat(n,1)=X(J);
end

end