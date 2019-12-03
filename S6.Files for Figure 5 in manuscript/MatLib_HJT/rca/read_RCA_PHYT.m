function [PHYT,StopInfo] = read_RCA_PHYT(FileToRead,GridInfo,iStep)
% function read_ecom_ts a function to read scalar variable obtained from
% ECOM simulations
% Created by Jiatang, SYSU
PHYT.time = [];
PHYT.PHYT1 = [];
PHYT.PHYT2 = [];
PHYT.PHYT3 = [];
% for iStep = iniStep:endStep
try
    m = memmapfile(...
        FileToRead,   ...
        'Format', {   ...
        'uint32'  [1]  'tag1';
        'single'  [1]  'time';
        'uint32'  [1]  'tag2';
        'uint32'  [1]  'tag3';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz-1] 'PHYT1';
        'uint32'  [1]  'tag4';
        'uint32'  [1]  'tag5';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz-1] 'PHYT2';
        'uint32'  [1]  'tag6';
        'uint32'  [1]  'tag7';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz-1] 'PHYT3';
        'uint32'  [1]  'tag8';
        }, ...
        'Repeat',1,  ...
        'offset',(iStep-1)*4*(3+3*(2+GridInfo.nx*GridInfo.ny*(GridInfo.nz-1))) ... %数据块的大小
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

PHYT.time = double(mdata(1).time);
PHYT.PHYT1 = double(mdata(1).PHYT1(:,:,:));
PHYT.PHYT2 = double(mdata(1).PHYT2(:,:,:));
PHYT.PHYT3 = double(mdata(1).PHYT3(:,:,:));
% end

return
end