%% main program to draw model results from ECOM
% Created by wangbin, SYSU
% SSS / SBS / SCHL / DO Bottom / DIC FLUX / PP / ASF / SOD
%
clear
clc
close all
%% Input model settings

% add matlib to ML's path
MatLib = 'F:\LB\碳循环\RCA_SFM_EOF\new\MatLib_HJT';
addpath(genpath(MatLib));


MatfileDir = 'matfile';
if ~exist(MatfileDir,'dir')
    mkdir(MatfileDir);
end

% grid dimension
GridInfo.nz = 17;  % sigma level
GridInfo.nx = 183; % 列
GridInfo.ny = 186; % 行


% load coastline（珠三角岸线文件）
CoastalLine = load('precoast.dat');
Coastal.lon = CoastalLine(:,1);
Coastal.lat = CoastalLine(:,2);
clearvars CoastalLine;


% Load grid data:
gridfile = fullfile('pregrid183by186.out');
griddata = load(gridfile);

GridInfo.lon = nan(GridInfo.nx,GridInfo.ny);
GridInfo.lat = nan(GridInfo.nx,GridInfo.ny);
GridInfo.depth = nan(GridInfo.nx,GridInfo.ny);
GridInfo.mask = nan(GridInfo.nx,GridInfo.ny);
GridInfo.angle = nan(GridInfo.nx,GridInfo.ny);
for irow = 1:length(griddata)
    jj = griddata(irow,1);  % 列
    ii = griddata(irow,2);  % 行
    GridInfo.lon(jj,ii) = griddata(irow,8);
    GridInfo.lat(jj,ii) = griddata(irow,7);
    GridInfo.depth(jj,ii) = griddata(irow,5);
    GridInfo.angle(jj,ii) = griddata(irow,6)*pi/180;
    if GridInfo.depth(jj,ii) < 0 || GridInfo.depth(jj,ii) > 9999
        fprintf('%4d%4d%8.1f \n',jj,ii,GridInfo.depth(jj,ii));
        error('Too small or too large values found in depth!');
    else
        GridInfo.mask(jj,ii) = 1;
    end
end
clearvars griddata;

%% Vertical Segmentation - Sigma Levels (KB)
SigmaLevel = [0.00  -0.01 -0.02  -0.05  -0.10  -0.20  -0.30  -0.40  -0.50 ...
    -0.60  -0.70  -0.80  -0.90  -0.95  -0.98  -0.99  -1.00];
DZZ = nan(1,numel(SigmaLevel)-1);
for ii = 1:numel(DZZ)
    DZZ(ii) = -(SigmaLevel(ii)+SigmaLevel(ii+1))/2;
end

Sigma_z = -diff(SigmaLevel);
delta_z = zeros(GridInfo.nx,GridInfo.ny,GridInfo.nz);
for kk = 1:GridInfo.nz-1
    delta_z(:,:,kk) = Sigma_z(kk)*GridInfo.depth;
end

grid_z =  zeros(GridInfo.nx,GridInfo.ny,GridInfo.nz);
for kk = 1:GridInfo.nz-1
    grid_z(:,:,kk) = DZZ(kk)*GridInfo.depth;
end

%% Cross-sections setting
SecMatFileName = 'Matfile\SectionGridPoints.mat';
if 1 %~exist(SecMatFileName,'file')
    % Section one
    SecLoc{1}.SecId = 'Transect A';
    SecLoc{1}.Grid_IX = 54*ones(numel([78:-1:35]),1);
    SecLoc{1}.Grid_IY = [78:-1:35];
    SecLoc{1}.NumPoint = numel(SecLoc{1}.Grid_IX);
    for ii = 1:SecLoc{1}.NumPoint
        ix = SecLoc{1}.Grid_IX(ii);
        iy = SecLoc{1}.Grid_IY(ii);
        SecLoc{1}.Lon(ii) = GridInfo.lon(ix,iy);
        SecLoc{1}.Lat(ii) = GridInfo.lat(ix,iy);
        SecLoc{1}.Bdep(ii) = GridInfo.depth(ix,iy);
        if ii==1
            SecLoc{1}.dis(ii) =  distance(SecLoc{1}.Lon(1),SecLoc{1}.Lat(1),SecLoc{1}.Lon(1),SecLoc{1}.Lat(1),6378.1);
        else
            SecLoc{1}.dis(ii) =  SecLoc{1}.dis(ii-1)+distance(SecLoc{1}.Lon(ii-1),SecLoc{1}.Lat(ii-1),SecLoc{1}.Lon(ii),SecLoc{1}.Lat(ii),6378.1);
        end
    end
    
    % Section two
    SecLoc{2}.SecId = 'Transect B';
    %     SecLoc{2}.Grid_IX = [51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51];
    %     SecLoc{2}.Grid_IY = [61 60 59 58 57 56 55 54 53 52 51 50 49 48 47 46 45];
    SecLoc{2}.Grid_IX = 107*ones(numel([170:-1:35]),1);
    SecLoc{2}.Grid_IY = [170:-1:35];
    SecLoc{2}.NumPoint = numel(SecLoc{2}.Grid_IX);  %=170个
    for ii = 1:SecLoc{2}.NumPoint
        ix = SecLoc{2}.Grid_IX(ii);
        iy = SecLoc{2}.Grid_IY(ii);
        SecLoc{2}.Lon(ii) = GridInfo.lon(ix,iy);
        SecLoc{2}.Lat(ii) = GridInfo.lat(ix,iy);
        SecLoc{2}.Bdep(ii) = GridInfo.depth(ix,iy);
        if ii==1
            SecLoc{2}.dis(ii) =  distance(SecLoc{2}.Lon(1),SecLoc{2}.Lat(1),SecLoc{2}.Lon(1),SecLoc{2}.Lat(1),6378.1);
        else
            SecLoc{2}.dis(ii) =  SecLoc{2}.dis(ii-1)+distance(SecLoc{2}.Lon(ii-1),SecLoc{2}.Lat(ii-1),SecLoc{2}.Lon(ii),SecLoc{2}.Lat(ii),6378.1);
        end
    end
    
    % Section three
    SecLoc{3}.SecId = 'Transect C';
    SecLoc{3}.Grid_IX = [66:133];
    SecLoc{3}.Grid_IY = 70*ones(numel([60:141]),1);
    SecLoc{3}.NumPoint = numel(SecLoc{3}.Grid_IX);
    for ii = 1:SecLoc{3}.NumPoint
        ix = SecLoc{3}.Grid_IX(ii);
        iy = SecLoc{3}.Grid_IY(ii);
        SecLoc{3}.Lon(ii) = GridInfo.lon(ix,iy);
        SecLoc{3}.Lat(ii) = GridInfo.lat(ix,iy);
        SecLoc{3}.Bdep(ii) = GridInfo.depth(ix,iy);
        if ii==1
            SecLoc{3}.dis(ii) =  distance(SecLoc{3}.Lon(1),SecLoc{3}.Lat(1),SecLoc{3}.Lon(1),SecLoc{3}.Lat(1),6378.1);
        else
            SecLoc{3}.dis(ii) =  SecLoc{3}.dis(ii-1)+distance(SecLoc{3}.Lon(ii-1),SecLoc{3}.Lat(ii-1),SecLoc{3}.Lon(ii),SecLoc{3}.Lat(ii),6378.1);
        end
    end
    
    % Section four
    save(SecMatFileName,'SecLoc');
else
    load(SecMatFileName);
end
%load('matfile/mycolor.mat');
%load('matfile/mycolor4.mat');
%% reading setting
% Time window for analysis
AnalysisModelTime{1} = [datenum(2006,07,01,0,0,0) datenum(2006,08,01,0,0,0)];   % July
%
AnalysisModelTime{2} = [datenum(2006,07,01,0,0,0) datenum(2006,08,01,0,0,0)];   % July to Augest
%
for month = [2];
    switch month
        case 12
            imonth = 'Dec';
            str = 'Dec';
        case 1
            imonth = 'Jan';
            str = 'Jan';
        case 2
            imonth = 'Feb';
            str = 'Feb';
        case 3
            imonth = 'Mar';
            str = 'Mar';
        case 4
            imonth = 'Apr';
            str = 'Apr';
        case 5
            imonth = 'May';
            str = 'May';
        case 6
            imonth = 'Jun';
            str = 'Jun';
        case 7
            imonth = 'Jul';
            str = 'Jul';
        case 8
            imonth = 'Aug';
            str = 'Aug';
        case 9
            imonth = 'Sep';
            str = 'Sep';
        case 10
            imonth = 'Oct';
            str = 'Oct';
        case 11
            imonth = 'Nov';
            str = 'Nov';
        case 13
            imonth = 'Summer';
            str = 'Summer';
        case 14
            imonth = 'Winter';
            str = 'Winter';
            %     case 1
            %         imonth = 'Dec';
            %         str = 'Winter';
            %     case 2
            %         imonth = 'Jan';
            %         str = 'Winter';
            %     case 3
            %         imonth = 'Feb';
            %         str = 'Winter';
            %     case 4
            %         imonth = 'Mar';
            %         str = 'Spring';
            %     case 5
            %         imonth = 'Apr';
            %         str = 'Spring';
            %     case 6
            %         imonth = 'May';
            %         str = 'Spring';
            %     case 7
            %         imonth = 'Jun';
            %         str = 'Summer';
            %     case 8
            %         imonth = 'Jul';
            %         str = 'Summer';
            %     case 9
            %         imonth = 'Aug';
            %         str = 'Summer';
            %     case 10
            %         imonth = 'Sep';
            %         str = 'Autumn';
            %     case 11
            %         imonth = 'Oct';
            %         str = 'Autumn';
            %     case 12
            %         imonth = 'Nov';
            %         str = 'Autumn';
    end
    CaseToShow = cell(1,numel(AnalysisModelTime));
    CaseToShow{1} = imonth;
    CaseToShow{2} = 'July to Augest';
    
    CasesSuffixAll = {'1999YES','1999NO'};
    
    CasesNamesAll = {'YES99',['2006 ',str]};
    
    
    % CasesSuffixAll = {'ERAsynopmonthly','ERAsynopdaily','ERAsynop6h','ERAhomo6h','ERAsynopVertical25_1-3D_m2',...
    %                'ERAsynopVertical25_1-3D_day',...
    %                'ERAsynopVertical25_3D_day',...
    %                'ERAsynopVertical25_1-3D_month',...
    %                'ERAsynopVertical25_3D_month'};
    % CasesNamesAll = {'#W2','#W1','Control','#W3','#T1','#D1','#D2','#D3','#D4'};
    SelectPlotId = 2;
    SelectPlotPeriods= [1];
    
    % for icase = 1:numel(CasesSuffixAll)
    %     CasesSuffix = CasesSuffixAll(icase);
    %     CasesName   = CasesNamesAll(icase);
    %     MatFileName = sprintf('JulyToAugest2006_%s.mat',CasesName{1});
    %     MatFileName = fullfile(MatfileDir,MatFileName);
    %     if ~exist(MatFileName,'file')
    %         fprintf('   -- reading %s \n', CasesName{1});
    %         eval('SubCode_Read');
    %         save(MatFileName,'Ecom3D','Coastal','GridInfo');
    %     end
    % end
    
    
    CasesSuffix = CasesSuffixAll(SelectPlotId);
    CasesName   = CasesNamesAll(SelectPlotId);
    MatFileName = sprintf('JulyToAugest2006_%s.mat',CasesName{1});
    MatFileName = fullfile(MatfileDir,MatFileName);
    % load(MatFileName);
    
    
    %% position subplots (2*2) and do interpolation
    for ii = 1:numel(SelectPlotPeriods)
        SelectPlotPeriod = SelectPlotPeriods(ii);
        FigureDir{ii} = sprintf('CombFigs_%s/%s',CasesName{1},CaseToShow{SelectPlotPeriod});
        if ~exist(FigureDir{ii},'dir')
            mkdir(FigureDir{ii});
        end
    end
    
    for ii = 1:numel(SelectPlotPeriods)
        SelectPlotPeriod = SelectPlotPeriods(ii);
        FigureDir{ii} = sprintf('CombFigs_%s',CasesName{1});
        if ~exist(FigureDir{ii},'dir')
            mkdir(FigureDir{ii});
        end
    end
    
    subPos = cell(1,2);
    % Positioning of all the subplots (2 by 2) for horizontal
    upM = .1;
    upMidM = .0;
    loMidM = .05;
    loM = .05;
    leftM = .05;
    rightMidM = .0;
    leftMidM = .05;
    rightM = .05;
    wC = (1 - leftM - rightM - rightMidM - leftMidM)/2;
    hC = (1 - upM - loM - upMidM - loMidM)/2;
    %topRow = loM + loMidM + upMidM + 2*hC;
    midRow = loM + loMidM + hC;
    botRow = loM;
    leftCol  = leftM;
    midCol   = leftM + leftMidM + wC;
    %rightCol = leftM + rightMidM + 2*wC;
    
    subPos{1} = [leftCol midRow wC hC;    midCol midRow wC hC;
        leftCol botRow wC hC;    midCol botRow wC hC;];
    
    % Positioning of all the subplots (3 by 2) for vertical
    upM = .1;
    upMidM = .1;
    loMidM = .1;
    loM = .1;
    leftM = .06;
    leftMidM = .05;
    rightMidM = .05;
    rightM = .02;
    wC = (1 - leftM - rightM - rightMidM)/2;
    hC = (1 - upM - loM - upMidM - loMidM)/3;
    topRow = loM + loMidM + upMidM + 2*hC;
    midRow = loM + loMidM + hC;
    botRow = loM;
    leftCol  = leftM;
    midCol   = leftM + leftMidM + wC;
    rightCol = leftM + rightMidM + 2*wC;
    
    % midCol   = leftM + wC*0.7;
    subPos{2} = [leftCol topRow wC hC;    midCol topRow wC hC;
        leftCol midRow wC hC;    midCol midRow wC hC;
        leftCol botRow wC hC;    midCol botRow wC hC;];
    
    %
    % To do interpolation
    %
    [NewLon,NewLat] = meshgrid(112:0.004:116,20:0.004:23);
    
    WetGrid.Index = find(~isnan(GridInfo.lon));
    WetGrid.Lon = GridInfo.lon(WetGrid.Index);
    WetGrid.Lat = GridInfo.lat(WetGrid.Index);
    
    %% Figure 1: Horizontal for Global
    IsToPlotHorizontal_1 = 0;
    
    PropertyOpt.xlim = [112.4,115.6];
    PropertyOpt.ylim = [20.9,22.8];
    
    if IsToPlotHorizontal_1
        varnames = {'PHYT','POC','DOC','pCO2'};
        varnames = {'pCO2'};
        for ivar = 1:numel(varnames)
            figure(ivar)
            switch varnames{ivar}
                case 'POC'
                    PropertyOpt.caxis{1} = [0 2.5];
                    PropertyOpt.xtick{1} = 0:0.5:2.5;
                    PropertyOpt.Colors{1} = 0:0.01:2.5;
                    PropertyOpt.ColorLines{1} = [0.1:0.1:3.6];
                    
                    PropertyOpt.caxis{2} = [0 2.5];
                    PropertyOpt.xtick{2} = 0:0.5:2.5;
                    PropertyOpt.Colors{2} = 0:0.01:2.5;
                    PropertyOpt.ColorLines{2} = [0.1:0.1:3.6];
                    
                    unit = 'POC (mg/L)';
                    load(['matfile/',imonth,'_POC.mat'])
                case 'DIC'
                    PropertyOpt.caxis{1} = [1.6 2.0];
                    PropertyOpt.xtick{1} = 1.6:0.1:2.0;
                    PropertyOpt.Colors{1} = 1.6:0.01:2.0;
                    PropertyOpt.ColorLines{1} = [1.6:0.05:2.0];
                    
                    PropertyOpt.caxis{2} = [1.6 2.0];
                    PropertyOpt.xtick{2} = 1.6:0.1:2.0;
                    PropertyOpt.Colors{2} = 1.6:0.01:2.0;
                    PropertyOpt.ColorLines{2} = [1.6:0.05:2.0];
                    unit = 'DIC (mmol/L)';
                    load(['matfile/',imonth,'_DIC.mat'])
                case 'DOC'
                    PropertyOpt.caxis{1} = [0 2.5];
                    PropertyOpt.xtick{1} = 0:0.5:2.5;
                    PropertyOpt.Colors{1} = 0:0.01:2.5;
                    PropertyOpt.ColorLines{1} = [0.2:0.4:3.6];
                    
                    PropertyOpt.caxis{2} = [0 2.5];
                    PropertyOpt.xtick{2} = 0:0.5:2.5;
                    PropertyOpt.Colors{2} = 0:0.01:2.5;
                    PropertyOpt.ColorLines{2} = [0.2:0.4:3.6];
                    unit = 'DOC (mg/L)';
                    load(['matfile/',imonth,'_DOC.mat'])
                case 'PHYT'
                    PropertyOpt.caxis{1} = [0 2.5]/10;
                    PropertyOpt.xtick{1} = [0:0.5:2.5]/10;
                    PropertyOpt.Colors{1} = [0:0.01:2.5]/10;
                    PropertyOpt.ColorLines{1} = [0.2:0.4:3.6]/10;
                    
                    PropertyOpt.caxis{2} = [0 2.5]/10;
                    PropertyOpt.xtick{2} = [0:0.5:2.5]/10;
                    PropertyOpt.Colors{2} = [0:0.01:2.5]/10;
                    PropertyOpt.ColorLines{2} = [0.2:0.4:3.6]/10;
                    unit = 'PHYT (mg/L)';
                    load(['matfile/',imonth,'_PHYT.mat'])
                case 'pCO2'
                    PropertyOpt.caxis{1} = [300 1500];
                    PropertyOpt.xtick{1} = 300:200:1500;
                    PropertyOpt.Colors{1} = 200:100:2000;
                    PropertyOpt.ColorLines{1} =[300:100:1500];   % [300:100:1500]
                    
                    PropertyOpt.caxis{2} = [300 1500];
                    PropertyOpt.xtick{2} = 300:200:1500;
                    PropertyOpt.Colors{2} = 200:100:2000;
                    PropertyOpt.ColorLines{2} = [300:100:1500];
                    unit = 'pCO2 (uatm)';
                    load(['matfile/',imonth,'_pCO2.mat'])
                otherwise
                    error('no such varname available: %s\n',varnames{ivar});
            end
            
            PropertyOpt.ColorbarPos{1} = [0.255 0.53 0.2 0.01];
            PropertyOpt.ColorbarPos{2} = [0.73 0.53 0.2 0.01];
            PropertyOpt.ColorbarPos{3} = [0.255 0.08 0.2 0.01];
            PropertyOpt.ColorbarPos{4} = [0.73 0.08 0.2 0.01];
            
            for iPeriod = 1:numel(SelectPlotPeriods)
                SelectPlotPeriod = SelectPlotPeriods(iPeriod);
                LayerToShow = {'surface layer','bottom layer','surface layer (std)','bottom layer (std)'};
                for isub = 1:2
                    h_sub = subplot(2,2,isub);
                    switch isub
                        case 1
                            eval(['OriginalModelData = ',varnames{ivar},'(:,:,1);']);
                            ind =1;
                        case 2
                            eval(['OriginalModelData = ',varnames{ivar},'(:,:,16);']);
                            ind =1;
                        case 3
                            str = sprintf('OriginalModelData = Ecom3D(SelectPlotPeriod).std.%s(:,:,1);',varnames{ivar});
                            eval(['OriginalModelData = ',varnames{ivar},'(:,:,1);']);
                            eval(str);
                            ind =2;
                        case 4
                            eval(['OriginalModelData = ',varnames{ivar},'(:,:,16);']);
                            ind =2;
                    end
                    %
                    % Do the interpoation onto a regular grid
                    %
                    %                 OriginalModelData = Ecom3D(1).salt(:,:,1);
                    GoodData = OriginalModelData(WetGrid.Index);
                    F = TriScatteredInterp(WetGrid.Lon,WetGrid.Lat,GoodData);
                    InterpData = F(NewLon,NewLat);
                    contourf(NewLon,NewLat,InterpData,PropertyOpt.Colors{ind},'LineStyle','none'); hold on
                    % surf(NewLon,NewLat,InterpData);shading interp;view(0,90);hold on
                    FillLandMap(Coastal.lon,Coastal.lat);   % 画岸线及其陆地填充
                    if ind==1
                        [c, h] = contour(GridInfo.lon,GridInfo.lat,OriginalModelData,PropertyOpt.ColorLines{ind}); hold on
                        set(h,'LineColor',[ 1 1 1]*1)
                        clabel(c, h);
                    end
                    load('mycolor4.mat')
                    colormap(mycolor4);
                    caxis(PropertyOpt.caxis{ind});
                    xlim(PropertyOpt.xlim);
                    ylim(PropertyOpt.ylim);
                    H1 = text(112.5,22.6,CasesName);
                    H2 = text(112.5,22.4,LayerToShow{isub});
                    H2 = text(112.5,22.2,CaseToShow{SelectPlotPeriod});
                    % H3 = text(114.45,21.05,unit);
                    amerc;
                    set(h_sub,'position',subPos{1}(isub,:));
                    han = colorbar('Location','NorthOutside','Position',PropertyOpt.ColorbarPos{isub});
                    set(get(han,'Title'),'string',unit,'HorizontalAlignment','left','VerticalAlignment','middle');
                    set(han,'xtick',PropertyOpt.xtick{ind},'TickDir','out','Box','on');
                end
                
                str = sprintf('%s/%s_%s_%s_11',FigureDir{iPeriod},CaseToShow{SelectPlotPeriod},varnames{ivar},CasesName{1});
                print('-dpng','-r400',str);
                print('-depsc',str);
                %   close;
            end
        end
    end
    %% Figure 2: Horizontal for local
    IsToPlotHorizontal_2 = 1;
    
    PropertyOpt.xlim = [112.9,114.1];
    PropertyOpt.ylim = [21.8,22.8];
    
    load('PearlIndex.mat')
    if IsToPlotHorizontal_2
        % SSS / SBS / SCHL / DO Bottom / DIC FLUX / PP / ASF / SOD
        varnames = {'PHYT','POC','DOC','pCO2','DIC','TPOC'};
        varnames = {'pCO2'};
        for ivar = 1:numel(varnames)
            figure(ivar)
            switch varnames{ivar}
                case 'POC'
                    PropertyOpt.caxis{1} = [0 2.5];
                    PropertyOpt.xtick{1} = 0:0.5:2.5;
                    PropertyOpt.Colors{1} = 0:0.05:2.5;
                    PropertyOpt.ColorLines{1} = [0:0.4:3.6];
                    PropertyOpt.ContourLines{1} = [0:0.4:3.6];
                    PropertyOpt.caxis{2} = [0 2.5];
                    PropertyOpt.xtick{2} = 0:0.5:2.5;
                    PropertyOpt.Colors{2} = 0:0.05:2.5;
                    PropertyOpt.ColorLines{2} = [0:0.4:3.6];
                    
                    unit = 'non-living POC (mg/L)';
                    load(['matfile/',imonth,'_POC.mat'])
                    fig_s = {'(g)','(h)'};
                case 'TPOC'
                    PropertyOpt.caxis{1} = [0 2.5];
                    PropertyOpt.xtick{1} = 0:0.5:2.5;
                    PropertyOpt.Colors{1} = 0:0.05:2.5;
                    PropertyOpt.ColorLines{1} = [0:0.4:3.6];
                    PropertyOpt.ContourLines{1} = [0:0.4:3.6];
                    PropertyOpt.caxis{2} = [0 2.5];
                    PropertyOpt.xtick{2} = 0:0.5:2.5;
                    PropertyOpt.Colors{2} = 0:0.05:2.5;
                    PropertyOpt.ColorLines{2} = [0:0.4:3.6];
                    fig_s = {'(e)','(f)'};
                    unit = 'POC (mg/L)';
                    load(['matfile/',imonth,'_POC.mat'])
                    load(['matfile/',imonth,'_PHYT.mat'])
                    TPOC = POC+PHYT;
                case 'DIC'
                    PropertyOpt.caxis{1} = [1.1 2];
                    PropertyOpt.xtick{1} = 1.1:0.3:2;
                    PropertyOpt.Colors{1} = 1.1:0.01:2;
                    PropertyOpt.ColorLines{1} = [1.1:0.15:2];
                    PropertyOpt.ContourLines{1} = [1.1:0.15:2];
                    PropertyOpt.caxis{2} = [1.1 2];
                    PropertyOpt.xtick{2} = 1.1:0.3:2;
                    PropertyOpt.Colors{2} = 1.1:0.01:2;
                    PropertyOpt.ColorLines{2} = [1.1:0.15:2];
                    PropertyOpt.ColorLines{2} = [1.1:0.15:2];
                    unit = 'DIC (mmol/L)';
                    load(['matfile/',imonth,'_DIC.mat'])
                    fig_s = {'(a)','(b)'};
                case 'DOC'
                    PropertyOpt.caxis{1} = [0 2.5];
                    PropertyOpt.xtick{1} = 0:0.5:2.5;
                    PropertyOpt.Colors{1} = 0:0.05:2.5;
                    PropertyOpt.ColorLines{1} = [0:0.4:3.6];
                    PropertyOpt.ContourLines{1} = [0:0.4:3.6];
                    PropertyOpt.caxis{2} = [0 2.5];
                    PropertyOpt.xtick{2} = 0:0.5:2.5;
                    PropertyOpt.Colors{2} = 0:0.05:2.5;
                    PropertyOpt.ColorLines{2} = [0:0.4:3.6];
                    unit = 'DOC (mg/L)';
                    load(['matfile/',imonth,'_DOC.mat'])
                    fig_s = {'(c)','(d)'};
                case 'PHYT'
                    PropertyOpt.caxis{1} = [0 3]/10;
                    PropertyOpt.xtick{1} = [0:1:3]/10;
                    PropertyOpt.Colors{1} = [0:0.05:3]/10;
                    PropertyOpt.ColorLines{1} = [0:0.8:3.6]/10;
                    PropertyOpt.ContourLines{1} = [0:0.8:3.6]/10;
                    PropertyOpt.caxis{2} = [0 3]/10;
                    PropertyOpt.xtick{2} = [0:1:3]/10;
                    PropertyOpt.Colors{2} = [0:0.05:3]/10;
                    PropertyOpt.ColorLines{2} = [0:0.8:3.6]/10;
                    unit = 'PHYT (mg/L)';
                    load(['matfile/',imonth,'_PHYT.mat'])
                    fig_s = {'(a)','(b)'};
                case 'pCO2'
                    PropertyOpt.caxis{1} = [300 1500];
                    PropertyOpt.xtick{1} = 300:400:1500;
                    PropertyOpt.Colors{1} = 200:10:2000;
                    PropertyOpt.ColorLines{1} = [500:500:1500];
                    PropertyOpt.ContourLines{1} = [500:500:1500];
                    PropertyOpt.caxis{2} = [300 1500];
                    PropertyOpt.xtick{2} = 300:400:1500;
                    PropertyOpt.Colors{2} = 200:10:2000;
                    PropertyOpt.ColorLines{2} = [500:500:1500];
                    fig_s = {'(a)','(b)'};
                    unit = 'pCO2 (uatm)';
                    load(['matfile/',imonth,'_pCO2.mat'])
                case 'pH'
                    PropertyOpt.caxis{1} = [7.1 8.3];
                    PropertyOpt.xtick{1} = 7.2:0.3:8.2;
                    PropertyOpt.Colors{1} = 7.1:0.0125:8.3;
                    PropertyOpt.ColorLines{1} = [7.2:0.05:8.0];
                    PropertyOpt.ContourLines{1} = [7.2:0.1:8.0];
                    PropertyOpt.caxis{2} = [7.1 8.3];
                    PropertyOpt.xtick{2} = 7.2:0.3:8.2;
                    PropertyOpt.Colors{2} = 7.1:0.0125:8.3;
                    PropertyOpt.ColorLines{2} = [7.2:0.05:8.0];
                    unit = 'pH';
                    load(['matfile/',imonth,'_pH.mat'])
                    pH = -log10(pH);
                    fig_s = {'(a)','(b)'};
                case 'ASF'
                    PropertyOpt.caxis{1} = [-200 50];
                    PropertyOpt.xtick{1} = -200:50:50;
                    PropertyOpt.Colors{1} =-250:5:10;
                    PropertyOpt.ColorLines{1} = [-200:40:50];
                    PropertyOpt.ContourLines{1} = [-200:40:50];
                    PropertyOpt.caxis{2} = [-200 50];
                    PropertyOpt.xtick{2} = -200:50:50;
                    PropertyOpt.Colors{2} =-250:5:50;
                    PropertyOpt.ColorLines{2} = [-200:40:50];
                    PropertyOpt.ContourLines{2} = [-200:40:50];
                    unit = 'ASF (mmol/m2/day)';
                    load(['matfile/',imonth,'_ASF.mat'])
                    ASF(:,:,end) = ASF(:,:,1);
                    fig_s = {'(a)','(b)'};
                case 'DO'
                    PropertyOpt.caxis{1} = [2 8];
                    PropertyOpt.xtick{1} = 2:2:8;
                    PropertyOpt.Colors{1} = 0:0.1:8;
                    PropertyOpt.ColorLines{1} = [0:2:8];
                    PropertyOpt.ContourLines{1} = [0:2:8];
                    PropertyOpt.caxis{2} = [2 8];
                    PropertyOpt.xtick{2} = 2:2:8;
                    PropertyOpt.Colors{2} = 0:0.1:8;
                    PropertyOpt.ColorLines{2} = [0:2:8];
                    unit = 'DO (mg/L)';
                    load(['matfile/',imonth,'_DO.mat'])
                    fig_s = {'(a)','(b)'};
                case 'Salt'
                    PropertyOpt.caxis{1} = [0 35];
                    PropertyOpt.xtick{1} = 5:10:30;
                    PropertyOpt.Colors{1} = 0:0.5:35;
                    PropertyOpt.ColorLines{1} = [0:5:35];
                    PropertyOpt.ContourLines{1} = [0:5:35];
                    PropertyOpt.caxis{2} = [0 35];
                    PropertyOpt.xtick{2} = 5:10:30;
                    PropertyOpt.Colors{2} = 0:0.5:35;
                    PropertyOpt.ColorLines{2} = [0:5:35];
                    unit = 'Salinity (psu)';
                    load(['matfile/',imonth,'_Salt.mat'])
                case 'JDIC'
                    PropertyOpt.caxis{1} = [-10 150];
                    PropertyOpt.xtick{1} = 0:40:150;
                    PropertyOpt.Colors{1} =-10:5:150;
                    PropertyOpt.ColorLines{1} = [-10:40:150];
                    PropertyOpt.ContourLines{1} = [-10:40:150];
                    PropertyOpt.caxis{2} = [-10 150];
                    PropertyOpt.xtick{2} = 0:40:150;
                    PropertyOpt.Colors{2} =-10:5:150;
                    PropertyOpt.ColorLines{2} = [-10:40:150];
                    PropertyOpt.ContourLines{2} = [-10:40:150];
                    unit = 'JDIC';
                    load(['matfile/',imonth,'_JDIC.mat'])
                    JDIC = repmat(JDIC,[1,1,16]);
                case 'GPP'
                    PropertyOpt.caxis{1} = [0 15];
                    PropertyOpt.xtick{1} = 0:5:15;
                    PropertyOpt.Colors{1} = 0:0.5:20;
                    PropertyOpt.ColorLines{1} = [0:3:15];
                    PropertyOpt.ContourLines{1} = [0:3:15];
                    PropertyOpt.caxis{2} = [0 15];
                    PropertyOpt.xtick{2} = 0:0.5:15;
                    PropertyOpt.Colors{2} = 0:0.5:15;
                    PropertyOpt.ColorLines{2} = [0:3:15];
                    unit = 'GPP';
                    load(['matfile/',imonth,'_GPP.mat'])
                case 'TGPP'
                    PropertyOpt.caxis{1} = [0 3000];
                    PropertyOpt.xtick{1} = 0:500:3000;
                    PropertyOpt.Colors{1} = 0:50:3000;
                    PropertyOpt.ColorLines{1} = [0:300:300];
                    PropertyOpt.ContourLines{1} = [0:300:300];
                    PropertyOpt.caxis{2} = [0 300];
                    PropertyOpt.xtick{2} = 0:50:300;
                    PropertyOpt.Colors{2} = 0:50:300;
                    PropertyOpt.ColorLines{2} = [0:300:300];
                    unit = 'GPP';
                    load(['matfile/',imonth,'_GPP.mat'])
                    TGPP = sum(GPP,3)*1000;
                    TGPP(:,:,16) = TGPP;
                otherwise
                    error('no such varname available: %s\n',varnames{ivar});
            end
            eval([varnames{ivar},' = mean(',varnames{ivar},',4);']);
            PropertyOpt.ColorbarPos{1} = [0.112 0.815 0.15 0.01];
            PropertyOpt.ColorbarPos{2} = [0.588 0.815 0.15 0.01];
            PropertyOpt.ColorbarPos{3} = [0.112 0.365 0.15 0.01];
            PropertyOpt.ColorbarPos{4} = [0.588 0.365 0.15 0.01];
            
            for iPeriod = 1:numel(SelectPlotPeriods)
                SelectPlotPeriod = SelectPlotPeriods(iPeriod);
                LayerToShow = {'surface layer','bottom layer','surface layer (tidal)','bottom layer (tidal)'};
                for isub = 1:2
                    h_sub = subplot(2,2,isub);
                    switch isub
                        case 1
                            eval(['OriginalModelData = ',varnames{ivar},'(:,:,1);']);
                            ind = 1;
                        case 2
                            eval(['OriginalModelData = ',varnames{ivar},'(:,:,16);']);
                            ind = 1;
                        case 3
                            eval(['OriginalModelData = ',varnames{ivar},'(:,:,1);']);
                            ind = 2;
                        case 4
                            eval(['OriginalModelData = ',varnames{ivar},'(:,:,16);']);
                            ind = 2;
                    end
                    %
                    % Do the interpoation onto a regular grid
                    %
                    
                    if ~strcmp(varnames{ivar},'ASF') && ~strcmp(varnames{ivar},'JDIC')
                        id = find(OriginalModelData<0);
                        OriginalModelData(id)=0;
                    end
                    if strcmp(varnames{ivar},'JDIC')
                        OriginalModelData = 1000*OriginalModelData;
                    end
                    GoodData = OriginalModelData(WetGrid.Index);
                    F = TriScatteredInterp(WetGrid.Lon,WetGrid.Lat,GoodData);
                    InterpData = F(NewLon,NewLat);
                    Data = log10(InterpData);
                    contourf(NewLon,NewLat,InterpData,PropertyOpt.Colors{ind},'LineStyle','none'); hold on
                    % surf(GridInfo.lon,GridInfo.lat,OriginalModelData);shading interp;view(0,90);hold on
                    FillLandMap(Coastal.lon,Coastal.lat);   % 画岸线及其陆地填充
                    load('mycolor4.mat')
                    colormap(mycolor4);
                    caxis(PropertyOpt.caxis{ind});
                    xlim(PropertyOpt.xlim);
                    ylim(PropertyOpt.ylim);
                    H1 = text(112.95,22.53,LayerToShow{isub});
                    H2 = text(112.95,22.44,strrep(CasesName,'_','\_'));
                    H3 = text(112.94,22.73,unit);
                    H4 = text(113.1,22.37,fig_s{isub});
                    % H4 = text(112.95,22.35,CaseToShow{SelectPlotPeriod});
                    amerc;
                    if strcmp(varnames{ivar},'DIC') && isub == 1
                        hold on; plot(SecLoc{1}.Lon,SecLoc{1}.Lat,'w','LineWidth',1.5);
                        hold on; plot(SecLoc{2}.Lon,SecLoc{2}.Lat,'w','LineWidth',1.5);
                        hold on; plot(SecLoc{3}.Lon,SecLoc{3}.Lat,'w','LineWidth',1.5);
                        text(113.98,22.0,['\fontsize{14}'...
                            '{\color{white}B}']);
                        text(113.68,22.32,['\fontsize{14}'...
                            '{\color{white}C}']);
                        text(113.6,21.9,['\fontsize{14}'...
                            '{\color{white}A}']);
                    end
                    OriginalModelData(find(OriginalModelData == Inf)) = 0.0;
                    [c, h] = contour(GridInfo.lon,GridInfo.lat,OriginalModelData,PropertyOpt.ContourLines{1},'Fill','on'); hold on
                    set(h,'LineColor',[0.5 0.5 0.5]*0)
                    clabel(c, h);
                    
                    set(h_sub,'position',subPos{1}(isub,:));
                    han = colorbar('Location','NorthOutside','Position',PropertyOpt.ColorbarPos{isub});
                    % set(get(han,'Title'),'string',unit,'VerticalAlignment','middle','HorizontalAlignment','right');
                    set(han,'xtick',PropertyOpt.xtick{ind},'TickDir','out','Box','on');
                end
                
                str = sprintf('%s/%s_%s_%s_2',FigureDir{iPeriod},CaseToShow{SelectPlotPeriod},varnames{ivar},CasesName{1})
                print('-dpng','-r400',str);
                close;
            end
        end
    end
    
    %% Figure 3: Vertical for transections
    IsToPlotVertical = 0;
    
    if IsToPlotVertical
        varnames = {'DIC','PHYT','POC','DOC','pCO2'};
        varnames = {'DIC','POC','DOC','HYDSAL'};
        for ivar = 1:numel(varnames)
            figure(ivar)
            switch varnames{ivar}
                case 'POC'
                    PropertyOpt.caxis = [0 2.5];
                    PropertyOpt.xtick = 0:0.5:2.5;
                    PropertyOpt.Colors = 0:0.05:2.5;
                    unit = 'POC (mg/L)';
                    PropertyOpt.ContourLines{1} = [0:0.4:3.6];
                    load(['matfile/',imonth,'_POC.mat'])
                    fig_s = {'(g)','(h)','(i)'};
                case 'DIC'
                    PropertyOpt.caxis = [1.6 1.9];
                    PropertyOpt.xtick = 1.6:0.05:1.9;
                    PropertyOpt.Colors = 1.6:0.02:1.9;
                    unit = 'DIC (mmol/L)';
                    PropertyOpt.ContourLines{1} = [1.6:0.03:1.9];
                    load(['matfile/',imonth,'_DIC.mat'])
                    fig_s = {'(a)','(b)','(c)'};
                case 'DOC'
                    PropertyOpt.caxis = [0 2.5];
                    PropertyOpt.xtick = 0:0.5:2.5;
                    PropertyOpt.Colors = 0:0.05:2.5;
                    unit = 'DOC (mg/L)';
                    PropertyOpt.ContourLines{1} = [0:0.4:3.6];
                    load(['matfile/',imonth,'_DOC.mat'])
                    fig_s = {'(d)','(e)','(f)'};
                case 'PHYT'
                    PropertyOpt.caxis = [0 2.5]/10;
                    PropertyOpt.xtick = [0:0.5:2.5]/10;
                    PropertyOpt.Colors = [0:0.05:2.5]/10;
                    unit = 'PHYT (mg/L)';
                    PropertyOpt.ContourLines{1} = [0:0.4:3.6]/10;
                    load(['matfile/',imonth,'_PHYT.mat'])
                case 'pCO2'
                    PropertyOpt.caxis = [300 1500];
                    PropertyOpt.xtick = 300:400:1500;
                    PropertyOpt.Colors = 200:10:2000;
                    PropertyOpt.ColorLines{1} = [500:500:1500];
                    PropertyOpt.ContourLines{1} = [500:500:1500];
                    unit = 'pCO2 (uatm)';
                    load(['matfile/',imonth,'_pCO2.mat'])
                case 'HYDSAL'
                    PropertyOpt.caxis = [0 35];
                    PropertyOpt.xtick = 0:5:35;
                    PropertyOpt.Colors = 0:0.05:35;
                    PropertyOpt.ColorLines{1} = [0:5:35];
                    PropertyOpt.ContourLines{1} = [0:5:35];
                    unit = 'Salinity (psu)';
                    load(['matfile/',imonth,'_HYDSAL.mat'])
                    HYDSAL = mean(HYDSAL,4);
                    fig_s = {'(d)','(e)','(f)'};
                otherwise
                    error('no such varname available: %s\n',varnames{ivar});
            end
            
            PropertyOpt.ColorbarPos{1} = [0.077 0.715 0.21 0.01];
            PropertyOpt.ColorbarPos{2} = [0.562 0.715 0.21 0.01];
            PropertyOpt.ColorbarPos{3} = [0.077 0.415 0.21 0.01];
            PropertyOpt.ColorbarPos{4} = [0.562 0.415 0.21 0.01];
            PropertyOpt.ColorbarPos{5} = [0.077 0.115 0.21 0.01];
            PropertyOpt.ColorbarPos{6} = [0.562 0.115 0.21 0.01];
            PropertyOpt.XlabelPos{1}= [0.48 0.678 0.02 0.02];
            PropertyOpt.XlabelPos{2}= [NaN NaN NaN NaN];
            PropertyOpt.XlabelPos{3}= [0.48 0.378 0.02 0.02];
            PropertyOpt.XlabelPos{4}= [NaN NaN NaN NaN];
            PropertyOpt.XlabelPos{5}= [0.48 0.078 0.02 0.02];
            PropertyOpt.XlabelPos{6}= [NaN NaN NaN NaN];
            PropertyOpt.YlabelPos{1}= [0.02 0.90 0.02 0.02];
            PropertyOpt.YlabelPos{2}= [NaN NaN NaN NaN];
            PropertyOpt.YlabelPos{3}= [0.02 0.60 0.02 0.02];
            PropertyOpt.YlabelPos{4}= [NaN NaN NaN NaN];
            PropertyOpt.YlabelPos{5}= [0.02 0.30 0.02 0.02];
            PropertyOpt.YlabelPos{6}= [NaN NaN NaN NaN];
            PropertyOpt.UlabelPos{1}= [0.065 0.77 0.02 0.02];
            PropertyOpt.UlabelPos{2}= [NaN NaN NaN NaN];
            PropertyOpt.UlabelPos{3}= [0.065 0.47 0.02 0.02];
            PropertyOpt.UlabelPos{4}= [NaN NaN NaN NaN];
            PropertyOpt.UlabelPos{5}= [0.065 0.17 0.02 0.02];
            PropertyOpt.UlabelPos{6}= [NaN NaN NaN NaN];
            for iPeriod = 1:numel(SelectPlotPeriods)
                SelectPlotPeriod = SelectPlotPeriods(iPeriod);
                load(SecMatFileName);
                for isec = 1:3
                    h_sub = subplot(3,2,2*isec-1);
                    %                 str = sprintf('OriginalModelData = Ecom3D(SelectPlotPeriod).%s;',varnames{ivar});
                    %                 eval(str);
                    eval(['OriginalModelData = ',varnames{ivar},'(:,:,:);']);
                    OriginalModelData(:,:,17) = OriginalModelData(:,:,16);
                    %eval(['Regressioncoefficient.',varnames{ivar},'(:,:,:,isec) = OriginalModelData;']);
                    Value = [];
                    %                 load('E:\LB\RCA_SFM_EOF\RCAS.mat')
                    %                 if month == 3
                    %                     S = mean(S(:,:,:,31+31:31+31+28),4);
                    %                 elseif month == 8
                    %                     S = mean(S(:,:,:,31+31+28+31+30+31+30:31+31+28+31+30+31+30+31),4);
                    %                 end
                    for ipoint = 1:SecLoc{isec}.NumPoint
                        ix = SecLoc{isec}.Grid_IX(ipoint);
                        iy = SecLoc{isec}.Grid_IY(ipoint);
                        Value = [Value squeeze(OriginalModelData(ix,iy,1:end-1))];
                        eval(['Regressioncoefficient.',varnames{ivar},'_',num2str(isec),'(ipoint,:) = squeeze(OriginalModelData(ix,iy,1:end-1));']);
                    end
                    id = find(Value<0);
                    Value(id)=0;
                    Xmat = repmat(SecLoc{isec}.dis,GridInfo.nz-1,1);
                    Ymat = -DZZ'*SecLoc{isec}.Bdep;
                    %                 if isec == 3
                    %                     id = find(isnan(SecLoc{4}.Bdep));
                    %                     if ~isempty(id)
                    %                         SecLoc{4}.Bdep(id) = -9999;
                    %                         SecLoc{4}.dis(id) = linspace(26.7600,47.6000,24);
                    %                     end
                    %                 end
                    contourf(Xmat,Ymat,Value,PropertyOpt.Colors,'LineStyle','none'); hold  on;
                    %                 if isec==3
                    %                     surf(Xmat,Ymat,Value);shading interp;view(0,90);hold on
                    %                 end
                    [c, h] = contour(Xmat,Ymat,Value,PropertyOpt.ContourLines{1}); hold on
                    set(h,'LineColor',[1 1 1]*1)
                    clabel(c, h);
                    load('mycolor4.mat')
                    colormap(mycolor4);
                    caxis(PropertyOpt.caxis);
                    hold on;
                    fill([SecLoc{isec}.dis,fliplr(SecLoc{isec}.dis)],[min(-SecLoc{isec}.Bdep)*ones(size(SecLoc{isec}.Bdep)),fliplr(-SecLoc{isec}.Bdep)],1*[1 1 1],'EdgeColor','none')
                    plot(SecLoc{isec}.dis,-SecLoc{isec}.Bdep,'k','Linewidth',1)
                    xlim([0 max(SecLoc{isec}.dis)]);
                    ylim([min(-SecLoc{isec}.Bdep) 0]);
                    title(strcat(SecLoc{isec}.SecId))
                    set(h_sub,'position',subPos{2}(2*isec-1,:));
                    han = colorbar('Location','NorthOutside','Position',PropertyOpt.ColorbarPos{2*isec-1});
                    % set(get(han,'Title'),'string',unit,'VerticalAlignment','middle','HorizontalAlignment','right');
                    set(han,'xtick',PropertyOpt.xtick,'TickDir','out','Box','on');
                    if isec == 1
                        text(2,-15,fig_s{isec});
                    elseif isec == 2
                        text(4,-15,fig_s{isec});
                    elseif isec == 3
                        text(4,-12,fig_s{isec});
                    end
                    
                end
                for isec=1:3
                    annotation('textbox',PropertyOpt.UlabelPos{isec*2-1},'String',unit,'FitBoxToText','on','LineStyle','none');
                    annotation('textbox',PropertyOpt.XlabelPos{isec*2-1},'String',{'km'},'FitBoxToText','on','LineStyle','none');
                    annotation('textbox',PropertyOpt.YlabelPos{isec*2-1},'String',{'m'},'FitBoxToText','on','LineStyle','none');
                end
                picstr = sprintf('%s/%s_%s_%s_3',FigureDir{iPeriod},CaseToShow{SelectPlotPeriod},varnames{ivar},CasesName{1});
                print('-dpng','-r400',picstr);
                %             mail2me(CasesName{1},[CaseToShow{SelectPlotPeriod},varnames{ivar}],[picstr '.png'])
                %    close
            end
        end
    end
end