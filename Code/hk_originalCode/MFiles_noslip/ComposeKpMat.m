function KpMat=ComposeKpMat(X,N)

KpMat=zeros(2,N);

for n=1:N
    for i=1:2
        J=4*N+2*(n-1)+i;
        KpMat(i,n)=X(J);
    end
end

end