function [Limite,StopInfo] = read_RCA_Limite(FileToRead,GridInfo,iStep)
% function read_ecom_ts a function to read scalar variable obtained from
% ECOM simulations
% Created by Jiatang, SYSU
Limite.time = [];
Limite.limite1 = [];
Limite.limite2 = [];
Limite.limite3 = [];
% for iStep = iniStep:endStep
try
    m = memmapfile(...
        FileToRead,   ...
        'Format', {   ...
        'uint32'  [1]  'tag1';
        'single'  [1]  'time';
        'uint32'  [1]  'tag2';
        'uint32'  [1]  'tag3';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz-1] 'limite1';
        'uint32'  [1]  'tag4';
        'uint32'  [1]  'tag5';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz-1] 'limite2';
        'uint32'  [1]  'tag6';
        'uint32'  [1]  'tag7';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz-1] 'limite3';
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

Limite.time = double(mdata(1).time);
Limite.limite1 = double(mdata(1).limite1(:,:,:));
Limite.limite2 = double(mdata(1).limite2(:,:,:));
Limite.limite3 = double(mdata(1).limite3(:,:,:));
% end

return
end