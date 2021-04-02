function radiiVec=ComposeradiiVec_3(nlay,polythick)

radiiVec0=linspace(0,14e-3,nlay);
radiiVec1=zeros(fix(nlay/2),1);
radiiVec=zeros(fix(nlay/2),1);
count=0;
for i=1:nlay
    if rem(i,2)==0
    count=count+1;
    radiiVec1(count,1)=radiiVec0(1,i);
    end
end
radiiVec1=radiiVec1+polythick;
clear count


count=0;
for i=1:(nlay+fix(nlay/2))
    if rem(i,3)==0
        radiiVec(i,1)=radiiVec1(i/3,1);
    else
        count=count+1;
        radiiVec(i,1)=radiiVec0(1,count);
        
    end
end
clear count


end