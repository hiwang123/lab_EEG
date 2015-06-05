function cutPoint=findmin(data)
    len=length(data);
    M=mean(data);
    aaa=[];
    cutPoint=[];
    for i=2:len-1
        if(data(i)>M)
            if(length(aaa)~=0&min(data(aaa))<30)
               cutPoint=[cutPoint min(aaa)];
            end
            aaa=[];
        else if(data(i)<data(i-1)&data(i)<data(i+1)&data(i)<M)
            aaa=[aaa i];
        end
    end
end