function [Data]=preprocess_fm(rd, datamin, datamax)
    %RawData = load('Data2_Linear.dat');
    
    
    % Each row in data represents an object [x1, x2, x3, ..., xd, y], where x's
    % are features and y is the attribute to be predicted.
    
    RawData = rd; %[rd.X rd.Y];

    [DataRow, DataCol] = size(RawData);




    RawData_min = datamin; %min(RawData,[],1);
    RawData_max = datamax; %max(RawData,[],1);

    Data = (RawData-ones(DataRow,1)*RawData_min)./(ones(DataRow,1)*(RawData_max-RawData_min));  % Convert all attributes to [0, 1]
    Data = (Data-0.5).*2;                                                                       % Convert all attributes to [-1,1]

end