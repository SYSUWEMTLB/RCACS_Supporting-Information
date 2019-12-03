function [TS,StopInfo] = read_RCA_DO(FileToRead,GridInfo,iStep)
% function read_ecom_ts a function to read scalar variable obtained from
% ECOM simulations
% Created by Jiatang, SYSU

TS = [];
% for iStep = iniStep:endStep
try
    m = memmapfile(...
        FileToRead,   ...
        'Format', {   ...
        'uint32'  [1]  'tag1';
        'single'  [1]  'time';
        'uint32'  [1]  'tag2';
        'uint32'  [1]  'tag3';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz-1] 'TS';
        'uint32'  [1]  'tag4';
        }, ...
        'Repeat',1,  ...
        'offset',(iStep-1)*4*(3+2+GridInfo.nx*GridInfo.ny*(GridInfo.nz-1)) ... %数据块的大小
        );
    mdata = m.data;
catch ME
    StopInfo.StopStep = iStep-1;
    StopInfo.EndOfFile = 1;
    StopInfo.ME = ME;
    return;
end
StopInfo.StopStep = iStep;
StopInfo.EndOfFile = 0;

TS = mdata(1).TS(:,:,:);
%
TS = double(TS);
% end

return
end