function radiiVec=ComposeradiiVec_2(nlay,polythick,outerRadius)


alpha=zeros(nlay,1);
radiiVec=zeros((nlay*2)+2,1);
radiiVec_0=zeros(nlay+2,1);


for i=1:nlay
    if i==1
        alpha(i,1)=3/4;
    else
        alpha(i,1)=3/4+1/4*(alpha(i-1,1))^4;
    end
end

for i=1:(nlay+2)
    j=nlay+3-i;
    if j==nlay+2
        radiiVec_0(j,1)=outerRadius;
    elseif j==1
        radiiVec_0(j,1)=0;
    else
        radiiVec_0(j,1)=radiiVec_0(j+1,1)*alpha(j-1,1);
    end
end

for i=1:(2*nlay+2)
    if i==1
        radiiVec(i,1)=radiiVec_0(i,1);
    elseif rem(i,2)==0
        radiiVec(i,1)=radiiVec_0(i/2+1,1);
    else
        radiiVec(i,1)=radiiVec_0((i+1)/2,1)+polythick;
    end
end
        


end