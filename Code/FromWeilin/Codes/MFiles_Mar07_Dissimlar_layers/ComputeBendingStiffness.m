function EI=ComputeBendingStiffness(KMat, CMat,mMat,gMat,radiiVec,muMat,N)

EI=0;
for n=1:N
    an=radiiVec(n,1);
    bn=radiiVec(n+1,1);
    
    T1n=0;
    for i=1:4
        
        T1ni=KMat(i,n)*(CMat(1,3,n)+CMat(2,3,n)*(mMat(i,n)+1)-CMat(3,4,n)*gMat(i,n)*mMat(i,n));
        T1ni = T1ni*(an^(mMat(i,n)+2)-bn^(mMat(i,n)+2));
        T1ni=T1ni/(mMat(i,n)+2);    
    T1n=T1n+T1ni;
    end
     
    T2n=muMat(1,n)*(CMat(1,3,n)+3*CMat(2,3,n))-2*muMat(2,n)*CMat(3,4,n)-1;
    T2n=T2n*((an^4-bn^4)/4);
    
    
    EIn=(T1n+T2n);
    EIn=(pi/CMat(3,3,n))*EIn;
    
    EI=EI+EIn;
    
    clear  T1n T2n EIn;
end

end