
severe=[1,2,3,6,8,13,14,16,17,20,22,23,24,25];
moderate=[4,5,7,9,10,11,12,15,18,19,21,26];
for i =1:length(severe)
   load(strcat('E:\new2\',int2str(severe(i)),'\ratio.mat'));
   Ratio(i,:)=[ratio 1];
end

for i =1:length(moderate)
   load(strcat('E:\new2\',int2str(moderate(i)),'\ratio.mat'));
   Ratio2(i,:)=[ratio 0];
end

s_out=mean(Ratio);
m_out=mean(Ratio2);

save(strcat('../result/result.mat'),'Ratio','Ratio2');

plot(s_out(1:length(s_out)-1));
plot(m_out(1:length(m_out)-1));