RawData = load('Data2_Logistic.dat');
[DataRow, DataCol] = size (RawData);

% Each row in data represents an object [x1, x2, x3, ..., xd, y], where x's
% are features and y is the label in binary.





RawData_min = min(RawData,[],1);
RawData_max = max(RawData,[],1);

Data = (RawData-ones(DataRow,1)*RawData_min)./(ones(DataRow,1)*(RawData_max-RawData_min));  % Convert all attributes to [0, 1]
Data = [(Data(:,1:end-1)-0.5).*2, Data(:,end)];                                             % Convert X part to [-1, 1]; Y in binary {0, 1}



errSum = 0;

for rep = 1:500
    
    fold = rand(DataRow, 1);
    SepLine = (0<fold) & (fold<=0.2);
    Test = Data(SepLine,:);
    Train = Data(not(SepLine),:);
    
    % Data are seperated into Train part (80%) and Test part (20%)
   
    
    
    epsilon = 0.5;
    [w, b] = Functional_Logistic(Train, epsilon);
    
    % Call function
   
    
    
    errorRate = logClassification(Test, w, b)/size(Test, 1);
    errSum = errSum + errorRate;
    
    % Evaluation
end


disp(errSum/500);