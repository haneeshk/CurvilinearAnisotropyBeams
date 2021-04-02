function nuMat = ComposenuMat(X,N)

nuMat=zeros(1,N);
for n=1:N
    J=6*N+n;
    nuMat(1,n)=X(J);
end

end