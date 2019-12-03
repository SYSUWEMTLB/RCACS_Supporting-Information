% 定时执行程序
clear
clc

MatLib = 'D:\wangbin\MatLib_HJT';
addpath(genpath(MatLib))

% 
MainFolder = 'D:\wangbin\ECOM06\Hydro_cases\newTracer';
CaseNames = {'ERAsynopReaTracer_v2_negetive','ERAsynopReaTracer_v2_positive'};
SubFolder = '\Hydro\0607to0608\Link_prog_v2\';
FileName = 'startup';
isrun = ones(size(CaseNames));
ischeck = 1;
fprintf('The Matlab start to check %s\n',datestr(now))
while ischeck
    for icase = 1:numel(CaseNames)
        FileId = strcat(MainFolder,CaseNames{icase},SubFolder,FileName);
        if exist(FileId,'file')
            if isrun(icase) == 1
                str = sprintf('Dear Bin, your model case named %s has been done at %s.',CaseNames{icase},datestr(now));
                mail2me(CaseNames{icase},str);
                isrun(icase) = 0;
            end 
        end
    end
    pause(300) 
    fprintf('The Matlab is working %s\n',datestr(now))
    zeroId = find(isrun == 1);
    if isempty(zeroId)
        ischeck = 0;
    end
end
run('D:\wangbin\ECOM06\GIF\ModelAnalysisTracer.m')
run('D:\wangbin\ECOM06\waterage\calculateB.m')



