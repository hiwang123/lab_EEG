

load('../result/result.mat');
data=[Ratio; Ratio2];
data=data(randperm(size(data,1)),:);

dlen=size(data,1);
flen=size(data,2);

AC=[];
for i = 1:floor(dlen*0.2):dlen-floor(dlen*0.2)
    idxtest=[i:i+floor(dlen*0.2)-1];
    idxtrain=[];
    for j = 1:dlen
        if j<i || j>i+floor(dlen*0.2)-1
            idxtrain=[idxtrain j];
        end
    end
    

    %%train
    y=data(idxtrain,flen);
    x=data(idxtrain,1:flen-1);

    model= svmtrain(y, x ,['-s 0 -d 2 -t 1 -g 1 -r 1 -c 0.01']);
    [predicted_label, accuracy, decision_values] = svmpredict(y, x , model);


    % %%%test
    test_y=data(idxtest,flen);
    test_x=data(idxtest,1:flen-1);
    [predicted_label, accuracy, decision_values] = svmpredict(test_y, test_x, model);
    AC=[AC accuracy(1)/100];
end

result_AC=mean(AC);


