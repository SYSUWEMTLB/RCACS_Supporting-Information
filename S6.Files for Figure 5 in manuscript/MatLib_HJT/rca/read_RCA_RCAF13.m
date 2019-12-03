function [Data,StopInfo] = read_RCA_RCAF13(FileToRead,RCAInfo,iStep)
% 
% 
% 

Data = [];
NDMPS = RCAInfo.NDMPS;
NOSYS = RCAInfo.NOSYS;
TNums = RCAInfo.TNums;
DataNums = RCAInfo.DataNums;
% for iStep = iniStep:endStep
try
    m = memmapfile(...
        FileToRead ,...
        'Format',{...
        'single' [1 1] 'tag1';...
        'single' [5,NDMPS,NOSYS] 'Data';...
        'single' [(TNums-DataNums)/5,5] 'nouse';...
        'single' [1 1] 'tag2'},...
        'Writable', true, ...
        'Repeat',1, ...
        'offset',(iStep-1)*4*(1+TNums+1));
    mdata = m.data;
catch ME
    StopInfo.StopStep = iStep-1;
    StopInfo.EndOfFile = 1;
    StopInfo.ME = ME;
    return;
end
StopInfo.StopStep = iStep;
StopInfo.EndOfFile = 0;

Data = mdata(1).Data;
%
Data = double(Data);
% end

return
end