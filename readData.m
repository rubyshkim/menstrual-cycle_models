function [TIME,DATA,SEM] = readData(variableNames,datafile)
% variableNames: cell array containing names of variables (e.g., {'LH','FSH','E2','P4'})

N = length(variableNames);
DATA = cell(1,N); SEM = DATA;

for i = 1:N
    data = readcell(datafile,'Sheet',variableNames{i});
    TIME{i} = cell2mat(data(2:end,1));
    DATA{i} = cell2mat(data(2:end,2));
    SEM{i} = cell2mat(data(2:end,3));
end

end