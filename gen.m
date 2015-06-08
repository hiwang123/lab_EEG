
X=[];

for i=1:33
    x=load(strcat('../',int2str(i),'/ratio.mat'));
    X=[X; x.ratio(49:156)];
    %X=[X; x.ratio(49:84); x.ratio(121:156)];
end

X=X';
save('../result/allsbjt_ratio.mat','X');