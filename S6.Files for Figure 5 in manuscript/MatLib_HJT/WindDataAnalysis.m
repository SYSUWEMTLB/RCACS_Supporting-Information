% WindDataAnalysis
%********************************************************************************
% 
% analysis different wind products
%
% count mean speed, mean direction, and
% frequency of all direction wind for 
% wind products
%
% 注意：
%   selfCount方法没有剔除静风（即风速为0），但roseCount方法已将静风剔除，因此
% 在有静风的情况下，roseCount方法求出的平均风速更大。
%   使用时注意修改风场对照时间（RefTime=datenum('2005-11-01 00:00:0.0');）
%   注意修改要画风场的时间范围（CountInitialTime=datenum('2006-08-01 00:00:0.0'); CountEndTime=datenum('2006-08-31 23:00:0.0');）
%   注意修改self count下保存数据的excel文件名（基本就是windcount，可不改）
%   注意修改rose count下保存数据的excel文件名中的时间信息（基本是一个月一个文件，FileName='windCount_Augest'）
%   注意修改玫瑰图保存文件夹名字的时间信息（RoseFolder='Wind_Rose\Augest2006';）
%*********************************************************************************

close all
clear all
clc

%% Reading wind from wind file
% add matlib to ML's path
MatLib = 'D:\wangbin\ECOM\ECOM99\MatLib_HJT';
addpath(genpath(MatLib));

% reading wind from wind file
MainDir='D:\wangbin\ECOM\ECOM06\met_forcing\';
CaseSuffix={'NCEP_NCAR','HK_Airport','HK_Waglan','MC_Airport'};
CaseName  ={'NCEPEast','HKAirport','HKWaglan','MCAirport'};
MetforcingFile={'\NCEPEast_MetForcings_AprtoSep2006.out',...
                '\HKAirport_MetForcings_AprtoSep2006_Kd6.out'...
                '\Waglan_MetForcings_AprtoSep2006_Kd6.out',...
                '\MCAirport_MetForcings_AprtoSep2006_Kd6.out',...
                }; 

RefTime=datenum('2005-11-01 00:00:0.0'); % edit
% 
fprintf('Reading wind file...\n')
switchSelfCount = 0; % 两种统计风场的方法，统计结果在有静风条件下会有不同
switchRoseCount = 1;

for icase=1:numel(CaseSuffix)
    FileName=strcat(MainDir,CaseSuffix{icase},MetforcingFile{icase});
    if ~exist(FileName,'file')
        error('errors:no such file exsit!')
    else
        fid=fopen(FileName,'r');
        nline=0;
        while ~feof(fid)
            nline=nline+1;
            tline=fgetl(fid);
            if mod(nline,2)==0
               vecline=str2num(tline);
               wind.speed{icase}(nline/2,1)=vecline(1);
               wind.direction{icase}(nline/2,1)=vecline(2);
            else
                vecline=str2num(tline);
                wind.time{icase}((nline+1)/2,1)=vecline/24+RefTime;
            end
        end
        fclose(fid);
    end
end

%% Count the wind (self)
fprintf('Counting wind...\n')
CountInitialTime=datenum('2006-07-01 00:00:0.0'); % edit
CountEndTime=datenum('2006-08-31 23:00:0.0'); % edit

if switchSelfCount
    for icase=1:numel(CaseName)
        index=find(wind.time{icase}>=CountInitialTime & wind.time{icase}<=CountEndTime);
        Countwind.time{icase}=wind.time{icase}(index);
        Countwind.speed{icase}=wind.speed{icase}(index);
        Countwind.direction{icase}=wind.direction{icase}(index);
        Countwind.meanSpeed{icase}=mean(Countwind.speed{icase});
        Countwind.meanDirection{icase}=mean(Countwind.direction{icase});
        WindDircVar={'N','NE','E','SE','S','SW','W','NW'};
        for ivar=1:numel(WindDircVar)
            str=sprintf('Countwind.%sdirection{%d}=0;',WindDircVar{ivar},icase); % 统计风频率
            eval(str);
            str=sprintf('Countwind.%sspeed{%d}=0;',WindDircVar{ivar},icase); % 各频率风风速
            eval(str);
        end
        for ihour=1:numel(Countwind.time{icase})
            dirc=Countwind.direction{icase}(ihour);
            speed=Countwind.speed{icase}(ihour);
            % 用dircIndex来确定此时的风向,只包括下边界，不包括上边界。如：
            % 0°--> N
            % 22.5° --> NE
            dircIndex=floor((dirc+22.5)/45)+1;
            if dircIndex==9
                dircIndex=1;
            end
            % 统计各风向出现次数
            str=sprintf('Countwind.%sdirection{%d}=Countwind.%sdirection{%d}+1;',...
                WindDircVar{dircIndex},icase,WindDircVar{dircIndex},icase);
            eval(str);
            % 统计各风向风速
            str=sprintf('Countwind.%sspeed{%d}=Countwind.%sspeed{%d}+speed;',...
                WindDircVar{dircIndex},icase,WindDircVar{dircIndex},icase);
            eval(str);
            % 记录每时刻风向，以便与检验
            str=sprintf('Countwind.dircstr{%d}{%d}=''%s'';',...
                icase,ihour,WindDircVar{dircIndex});
            eval(str);
            
        end
        % 统计各风向频率和风速
        for ivar=1:numel(WindDircVar)
            str=sprintf('Countwind.%sfrequency{%d}=Countwind.%sdirection{%d}/%d;',...
                WindDircVar{ivar},icase,WindDircVar{ivar},icase,numel(Countwind.time{icase}));
            eval(str);
            str=sprintf('Countwind.%sspeed{%d}=Countwind.%sspeed{%d}/Countwind.%sdirection{%d};',...
                WindDircVar{ivar},icase,WindDircVar{ivar},icase,WindDircVar{ivar},icase);
            eval(str);
        end
    end
end

% Writing into a excel file
fprintf('Writing wind...\n')
FileName='windCount';
if switchSelfCount
    sheetName=sprintf('WindFrequency%sto%s',datestr(CountInitialTime,28),datestr(CountEndTime,28)); % edit
    Format='dd/mm/yyyy';
    TimeInfo=[datestr(CountInitialTime,Format) ' : ' datestr(CountEndTime,Format)];
    % time info
    xlswrite(FileName,{'Time Interval'},sheetName,'A1:A1');
    xlswrite(FileName,{TimeInfo},sheetName,'B1:B1');
    % case info
    xlsLabel={'B','C','D','E','F','G','H','I','J','K'};
    Range=sprintf('B2:%s2',xlsLabel{numel(CaseName)});
    xlswrite(FileName,CaseName,sheetName,Range);
    % label
    xlswrite(FileName,{'mean speed'},sheetName,'A3:A3');
    xlswrite(FileName,{'mean direction'},sheetName,'A4:A4');
    xlswrite(FileName,{'N WIND'},sheetName,'A5:A5');
    xlswrite(FileName,{'NE WIND'},sheetName,'A6:A6');
    xlswrite(FileName,{'E WIND'},sheetName,'A7:A7');
    xlswrite(FileName,{'SE WIND'},sheetName,'A8:A8');
    xlswrite(FileName,{'S WIND'},sheetName,'A9:A9');
    xlswrite(FileName,{'SW WIND'},sheetName,'A10:A10');
    xlswrite(FileName,{'W WIND'},sheetName,'A11:A11');
    xlswrite(FileName,{'NW WIND'},sheetName,'A12:A12');
    % count result
    for icase=1:numel(CaseName)
        % mean speed
        Range=sprintf('%s3:%s3',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.meanSpeed{icase},sheetName,Range);
        % mean direction
        Range=sprintf('%s4:%s4',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.meanDirection{icase},sheetName,Range);
        % N wind
        Range=sprintf('%s5:%s5',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.Nfrequency{icase},sheetName,Range);
        % NE wind
        Range=sprintf('%s6:%s6',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.NEfrequency{icase},sheetName,Range);
        % E wind
        Range=sprintf('%s7:%s7',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.Efrequency{icase},sheetName,Range);
        % SE wind
        Range=sprintf('%s8:%s8',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.SEfrequency{icase},sheetName,Range);
        % S wind
        Range=sprintf('%s9:%s9',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.Sfrequency{icase},sheetName,Range);
        % SW wind
        Range=sprintf('%s10:%s10',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.SWfrequency{icase},sheetName,Range);
        % W wind
        Range=sprintf('%s11:%s11',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.Wfrequency{icase},sheetName,Range);
        % NW wind
        Range=sprintf('%s12:%s12',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.NWfrequency{icase},sheetName,Range);
    end
    fprintf('%s/%s has been done! \n', FileName,sheetName);
    
    sheetName=sprintf('WindSpeed%sto%s',datestr(CountInitialTime,28),datestr(CountEndTime,28)); % edit
    Format='dd/mm/yyyy';
    TimeInfo=[datestr(CountInitialTime,Format) ' : ' datestr(CountEndTime,Format)];
    % time info
    xlswrite(FileName,{'Time Interval'},sheetName,'A1:A1');
    xlswrite(FileName,{TimeInfo},sheetName,'B1:B1');
    % case info
    xlsLabel={'B','C','D','E','F','G','H','I','J','K'};
    Range=sprintf('B2:%s2',xlsLabel{numel(CaseName)});
    xlswrite(FileName,CaseName,sheetName,Range);
    % label
    xlswrite(FileName,{'mean speed'},sheetName,'A3:A3');
    xlswrite(FileName,{'mean direction'},sheetName,'A4:A4');
    xlswrite(FileName,{'N WIND'},sheetName,'A5:A5');
    xlswrite(FileName,{'NE WIND'},sheetName,'A6:A6');
    xlswrite(FileName,{'E WIND'},sheetName,'A7:A7');
    xlswrite(FileName,{'SE WIND'},sheetName,'A8:A8');
    xlswrite(FileName,{'S WIND'},sheetName,'A9:A9');
    xlswrite(FileName,{'SW WIND'},sheetName,'A10:A10');
    xlswrite(FileName,{'W WIND'},sheetName,'A11:A11');
    xlswrite(FileName,{'NW WIND'},sheetName,'A12:A12');
    % count result
    for icase=1:numel(CaseName)
        % mean speed
        Range=sprintf('%s3:%s3',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.meanSpeed{icase},sheetName,Range);
        % mean direction
        Range=sprintf('%s4:%s4',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.meanDirection{icase},sheetName,Range);
        % N wind
        Range=sprintf('%s5:%s5',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.Nspeed{icase},sheetName,Range);
        % NE wind
        Range=sprintf('%s6:%s6',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.NEspeed{icase},sheetName,Range);
        % E wind
        Range=sprintf('%s7:%s7',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.Espeed{icase},sheetName,Range);
        % SE wind
        Range=sprintf('%s8:%s8',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.SEspeed{icase},sheetName,Range);
        % S wind
        Range=sprintf('%s9:%s9',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.Sspeed{icase},sheetName,Range);
        % SW wind
        Range=sprintf('%s10:%s10',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.SWspeed{icase},sheetName,Range);
        % W wind
        Range=sprintf('%s11:%s11',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.Wspeed{icase},sheetName,Range);
        % NW wind
        Range=sprintf('%s12:%s12',xlsLabel{icase},xlsLabel{icase});
        xlswrite(FileName,Countwind.NWspeed{icase},sheetName,Range);
    end
    fprintf('%s/%s has been done! \n', FileName,sheetName);
end

%% rosecount
if switchRoseCount
    RoseFolder='Wind_RoseTT\Augest2006';
    FileName='windCount';
    if ~exist(RoseFolder,'dir')
        mkdir(RoseFolder)
    end
    for icase=1:numel(CaseName)
        index=find(wind.time{icase}>=CountInitialTime & wind.time{icase}<=CountEndTime);
        TitleName=sprintf('%s in %s',CaseName{icase},datestr(CountEndTime,28));
        FigName=sprintf('%s from %s to %s',CaseName{icase},datestr(CountInitialTime,28),datestr(CountEndTime,28));
        [figure_handle,count,speeds,directions,Table] = ...
            WindRose_ORI(wind.direction{icase}(index),wind.speed{icase}(index),'titlestring',{TitleName;' '});
        FigName=fullfile(RoseFolder,FigName);
        print('-dpng','-r400',FigName)
        
%         % write into xlsx file
%         sheetName=sprintf('%s%sto%s',CaseName{icase},datestr(CountInitialTime,28),datestr(CountEndTime,28));
%         Format='dd/mm/yyyy';
%         TimeInfo=[datestr(CountInitialTime,Format) ' : ' datestr(CountEndTime,Format)];
%         % time info
%         xlswrite(FileName,{'Time Interval'},sheetName,'A1:A1');
%         xlswrite(FileName,{TimeInfo},sheetName,'B1:B1');
%         % case info
%         xlswrite(FileName,{'Case Name'},sheetName,'A2:A2');
%         xlswrite(FileName,{CaseName{icase}},sheetName,'B2:B2');
%         % wind info
%         [m,n]=size(Table);
%         xlsLabel={'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W'};
%         Range=sprintf('A3:%s%d',xlsLabel{n},m+2);
%         xlswrite(FileName,Table,sheetName,Range);
    end
end
fprintf('Job Done');

%%