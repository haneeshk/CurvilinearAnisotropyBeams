function KMat = ComposeKMat(X,N)

KMat=zeros(4,N);
for n=1:N
    for i=1:4
        J=4*(n-1)+i;
        KMat(i,n)=X(J);
    end
end

end
