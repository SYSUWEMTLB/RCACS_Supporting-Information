% 3D plot

clear all
close all
clc
MatLib = 'G:\LB\碳循环\RCA_SFM_EOF\new\MatLib_HJT';
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

if 1
    GridInfo.depth(73:75,44:45) = 2;
    GridInfo.depth(75,45) = 2;
    
    GridInfo.depth(69:71,46) = 2;
    
    GridInfo.depth(79:81,45:46) = 2;
    GridInfo.depth(82,46) = 2;
    
    GridInfo.depth(74:77,50) = 2;
    
    GridInfo.depth(78,56) = 2; %
    
    GridInfo.depth(71,52) = 2;
    
    % GridInfo.depth(63:67,62:64) = 2;
    GridInfo.depth(64:66,62:64) = 2;
    
    GridInfo.depth(91:93,46) = 2;
    
    GridInfo.depth(110:111,47) = 2;
    
    GridInfo.depth(113:115,35) = 2;
    
    GridInfo.depth(119,47) = 2;
    
    GridInfo.depth(122,45) = 2;
    
    GridInfo.depth(118:120,36) = 2;
    
    GridInfo.depth(127:128,48) = 2;
    
    
    % GridInfo.depth(55,59:72) = 2;
    % GridInfo.depth(56,59:66) = 2;
    % GridInfo.depth(57:60,60:66) = 2;
    % GridInfo.depth(59:60,59) = 2;
    % GridInfo.depth(61,61:66) = 2;
    GridInfo.depth(56:59,61:65) = 2;
    
    % GridInfo.depth(80:83,93:97) = 2;
    % GridInfo.depth(84:85,95:97) = 2;
    GridInfo.depth(81:82,94:96) = 2;
    GridInfo.depth(83:84,96) = 2;
    
    GridInfo.depth(51:52,62:71) = 2;
    
    GridInfo.depth(92,59) = 2;
    
    GridInfo.depth(97:98,59) = 2;
    
    GridInfo.depth(100:101,55:57) = 2;
    
    GridInfo.depth(110:111,53) = 2;
    
    GridInfo.depth(115:116,55) = 2;
    
    GridInfo.depth(117,57) = 2;
    
    GridInfo.depth(133,60) = 2;
    
    GridInfo.depth(128:129,58:59) = 2;
    
    GridInfo.depth(128,56) = 2;
    
    GridInfo.depth(131,56) = 2;
    
    GridInfo.depth(137,[52,55]) = 2;
    GridInfo.depth(138,52:55) = 2;
    
    GridInfo.depth(139:140,57:59) = 2;
    GridInfo.depth(140:141,55:58) = 2;
    GridInfo.depth(142,52:58) = 2;
    GridInfo.depth(143,51:57) = 2;
    GridInfo.depth(144,54:55) = 2;
    
    GridInfo.depth(139,66:68) = 2;
    
    GridInfo.depth(137,68) = 2;
    
    GridInfo.depth(39:40,51:54) = 2;
    GridInfo.depth(41,52:54) = 2;
    
    GridInfo.depth(13,[43,44,47:49]) = 2;
    GridInfo.depth(14,[40:44,48:50]) = 2;
    GridInfo.depth(15,[44:49]) = 2;
    GridInfo.depth(16,[45,47:51]) = 2;
    GridInfo.depth(17,49:51) = 2;
    GridInfo.depth(18,[40,50,51]) = 2;
    GridInfo.depth(19,50:51) = 2;
    
    GridInfo.depth([34,36],50) = 2;
    GridInfo.depth(34:36,51) = 2;
    
    GridInfo.depth(29,54:55) = 2;
    GridInfo.depth(35,56) = 2;
    
    GridInfo.depth(122,45) = 2;
    
    GridInfo.depth(110:111,47) = 2;
    
    GridInfo.depth(119,47) = 2;
    
    GridInfo.depth(110:111,53) = 2;
    
    GridInfo.depth(115:116,55) = 2;
    
    GridInfo.depth(117,57) = 2;
    
    GridInfo.depth(29,54:55) = 2;
    
    GridInfo.depth(6,46:48) = 2;
    GridInfo.depth(7,47:48) = 2;
    GridInfo.depth(8,46:50) = 2;
    GridInfo.depth(9,47:51) = 2;
    GridInfo.depth(10:11,50:51) = 2;
end

if 0
    PZ = zeros(size(GridInfo.depth));
    for i = 5:size(GridInfo.lon,1)-4
        for j = 5:size(GridInfo.lon,2)-4
            Z4(1) = double(isnan(GridInfo.depth(i,j)));
            Z4(2) = double(isnan(GridInfo.depth(i+1,j)));
            Z4(3) = double(isnan(GridInfo.depth(i+1,j+1)));
            Z4(4) = double(isnan(GridInfo.depth(i+1,j-1)));
            Z4(5) = double(isnan(GridInfo.depth(i-1,j)));
            Z4(6) = double(isnan(GridInfo.depth(i-1,j+1)));
            Z4(7) = double(isnan(GridInfo.depth(i-1,j-1)));
            Z4(8) = double(isnan(GridInfo.depth(i,j+1)));
            Z4(9) = double(isnan(GridInfo.depth(i,j-1)));
            if sum(Z4) == 9
                PZ(i,j) = 1;
            end
        end
    end
    iPZ = find(PZ == 1);
    for i = 1:numel(iPZ)
        if isnan(GridInfo.depth(iPZ(i)))
            GridInfo.depth(iPZ(i)) = 2;
        end
    end
end
GridInfo.lon2 = GridInfo.lon;
GridInfo.lat2 = GridInfo.lat;
GridInfo.lon(3:end-3,3:end-3) = interp_lb(GridInfo.lon(3:end-3,3:end-3));
GridInfo.lat(3:end-3,3:end-3) = interp_lb(GridInfo.lat(3:end-3,3:end-3));
GridInfo.depth(3:end-3,3:end-3) = interp_lb(GridInfo.depth(3:end-3,3:end-3));
GridInfo.lon1 = GridInfo.lon;
GridInfo.lat1 = GridInfo.lat;

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
[NewLon,NewLat] = meshgrid(112:0.004:116,20:0.004:23);
WetGrid.Index = find(~isnan(GridInfo.lon));
WetGrid.Lon = GridInfo.lon(WetGrid.Index);
WetGrid.Lat = GridInfo.lat(WetGrid.Index);

OriginalModelData = GridInfo.depth;
GoodData = OriginalModelData(WetGrid.Index);
F = TriScatteredInterp(WetGrid.Lon,WetGrid.Lat,GoodData);
InterpDepth = F(NewLon,NewLat);
clear OriginalModelData GoodData F


load('G:\LB\碳循环\RCA_SFM_EOF\new\matfile\Dec_pH.mat');
Data = -log10(pH);
Data(find(Data>10000)) = nan;
Data(find(Data <= 0)) = nan;

load('Compare_06Jul_1202.mat')
Obs = Compare.talk;
Observing = Obs(:,1);
Simulating = Obs(:,2);
Station_Lon = Obs(:,3);
Station_Lat = Obs(:,4);
Station_Level = Obs(:,5);
Observing(find(Observing >10000)) = nan;
Observing(find(Observing <= 0)) = nan;
Simulating(find(Simulating >10000)) = nan;
Simulating(find(Simulating <= 0)) = nan;


figure(1)

X = repmat(GridInfo.lon,[1 1 16]);
Y = repmat(GridInfo.lat,[1 1 16]);
Z = -1*grid_z;
V = Data;
Z3 = nan(size(GridInfo.lon,1),size(GridInfo.lon,2),1001);
Z3(:,:,end) = 100;
for ii = 1:size(GridInfo.lon,1)
    for jj = 1:size(GridInfo.lon,2)
        Z1 = Z(ii,jj,:);
        X1 = V(ii,jj,:);
        Z1(find(isnan(Z1))) = 0;
        X1(find(isnan(Z1))) = 0;
        Z1 = Z1(:);
        X1 = X1(:);
        Xnew(ii,jj,:) = nan(1,1,numel([-100:0.1:0]));
        if numel(find(Z1~=0)) ~= 0
            if numel(find(isnan(X1))) ~= 16
                Xnew1 = interp1(Z1(1:end-1),X1,[round(min(Z1(1:end-1))*10)/10:0.1:round(max((Z1(1:end-1)))*10)/10]);
                Xnew(ii,jj,[abs(round(max((Z1(1:end-1)))*10))+1:abs(round(min(Z1(1:end-1))*10))+1]) = flip(Xnew1);
                i1 = find(~isnan(Xnew1));
                Xnew(ii,jj,1:abs(round(max((Z1(1:end-1)))*10))+1) = Xnew1(i1(end))/10;
            end
        end
        
        if GridInfo.depth(ii,jj) > 0
            Z2 = [GridInfo.depth(ii,jj)];
            index1_Z2 = round(Z2(1)*10);
            Z3(ii,jj,[index1_Z2+1:1000+1]) = [-1*index1_Z2/10:-0.1:-100]/10;
        else
            Z3(ii,jj,1:end) = [0:-0.1:-100]/10;
        end
    end
end
Z31 = Z3;
Z31(~isnan(Z3)) = 1;
Z31(isnan(Z3)) = 0;
PZ = nan(183,186);
PZ1 = nan(183,186);

%%
clear V
OriginalModelData = Z3(:,:,1000);
GoodData = OriginalModelData(WetGrid.Index);
F = TriScatteredInterp(WetGrid.Lon,WetGrid.Lat,GoodData);
InterpDepth_200 = F(NewLon,NewLat);
clear OriginalModelData GoodData F
V = Xnew;
figure(1)
subplot(1,2,1)
axis([112.5 115.6 20.5 23 -10 1])
view(25,30)
set(gca,'DataAspectRatio',[1.2 1 12/2.5])
hold on
load('G:\LB\碳循环\RCA_SFM_EOF\new\mycolor4.mat')
colormap(mycolor4);
for ndepth = [ ]
    NewLon1 = NewLon;
    NewLat1 = NewLon;
    if 0
        OriginalModelData = Data(:,:,ndepth);
        OriginalModelData = interp_lb(OriginalModelData);
        GoodData = OriginalModelData(WetGrid.Index);
        F = TriScatteredInterp(WetGrid.Lon,WetGrid.Lat,GoodData);
        InterpData = F(NewLon,NewLat);
        i1 = find(GridInfo.lon1 > 115);
        i2 = find(GridInfo.lon1 < 103.8);
        i3 = find(GridInfo.lat1 > 24);
        i4 = find(GridInfo.lat1 < 1.5);
        i5 = find(isnan(InterpDepth_200));
        InterpData(i1) = nan;
        InterpData(i2) = nan;
        InterpData(i3) = nan;
        InterpData(i4) = nan;
        InterpData(i5) = nan;
        NewLon1(i1) = nan;
        NewLon1(i2) = nan;
        NewLon1(i3) = nan;
        NewLon1(i4) = nan;
        NewLon1(i5) = nan;
        NewLat1(i1) = nan;
        NewLat1(i2) = nan;
        NewLat1(i3) = nan;
        NewLat1(i4) = nan;
        NewLat1(i5) = nan;
        PropertyOpt.Colors{1} = 7.1:0.0125:8.3;
        %PropertyOpt.Colors{1} = 1.3:0.0125:2.1;
        i1 = find(GridInfo.lon1 > 115);
        i2 = find(GridInfo.lon1 < 103.8);
        i3 = find(GridInfo.lat1 > 24);
        i4 = find(GridInfo.lat1 < 1.5);
        i5 = find(isnan( Z3(:,:,1000)));
        GridInfo.lon(i1) = nan;
        GridInfo.lon(i2) = nan;
        GridInfo.lon(i3) = nan;
        GridInfo.lon(i4) = nan;
        GridInfo.lon(i5) = nan;
        GridInfo.lat(i1) = nan;
        GridInfo.lat(i2) = nan;
        GridInfo.lat(i3) = nan;
        GridInfo.lat(i4) = nan;
        GridInfo.lat(i5) = nan;
        OriginalModelData(i1) = nan;
        OriginalModelData(i2) = nan;
        OriginalModelData(i3) = nan;
        OriginalModelData(i4) = nan;
        OriginalModelData(i5) = nan;
        [h] = pcolor(GridInfo.lon,GridInfo.lat,OriginalModelData); hold on
        shading flat
        shading interp
        zdata = ones(size(get(h,'XData')));
        set(h,'ZData',-GridInfo.depth+2)
        if 0
            [C,h] = contour(GridInfo.lon,GridInfo.lat,OriginalModelData); hold on
            hh = get(h,'Children');
            hold on
            for i=1:numel(hh)
                XD = get(hh(i),'XData');
                YD = get(hh(i),'YData');
                for izd = 1:numel(XD)
                    m1 = (abs(XD(izd)-GridInfo.lon(:)));
                    m2 = (abs(YD(izd)-GridInfo.lat(:)));
                    izd1 = find(m1 <= 0.01);
                    izd2 = find(m2 <= 0.01);
                    izd3 = intersect(izd1,izd2);
                    if numel(izd3) >= 1
                        zdata(izd,1) = GridInfo.depth(izd3(1));
                    else
                        zdata(izd,1) = nan;
                    end
                end
                set(hh(i), 'ZData',-zdata)
            end
            set(h,'LineColor','w')
        end
    end
end


for nobs = 1:size(Obs,1)
    add1 = floor(Station_Level(nobs)/6)*4;
    %     scatter3(GridInfo.lon(Station_Lon(nobs),Station_Lat(nobs)), ...
    %         GridInfo.lat(Station_Lon(nobs),Station_Lat(nobs)), ...
    %         -Station_Level(nobs)-add1, ...
    %         10,Observing(nobs),'filled');
    if Station_Lon(nobs) > 0
        switch floor(Station_Level(nobs)/6)
            case 0
                scatter3(GridInfo.lon(Station_Lon(nobs),Station_Lat(nobs)), ...
                    GridInfo.lat(Station_Lon(nobs),Station_Lat(nobs)), ...
                    -add1-1, ...
                    50,Observing(nobs),'filled','^','MarkerEdgeColor','k');
                hold on
            case 1
                scatter3(GridInfo.lon(Station_Lon(nobs),Station_Lat(nobs)), ...
                    GridInfo.lat(Station_Lon(nobs),Station_Lat(nobs)), ...
                    -add1-1, ...
                    40,Observing(nobs),'filled','MarkerEdgeColor','k');
                hold on
            case 2
                scatter3(GridInfo.lon(Station_Lon(nobs),Station_Lat(nobs)), ...
                    GridInfo.lat(Station_Lon(nobs),Station_Lat(nobs)), ...
                    -add1-1, ...
                    50,Observing(nobs),'filled','v','MarkerEdgeColor','k');
                hold on
        end
    end
end
i1 = find(GridInfo.lon1 > 115);
i2 = find(GridInfo.lon1 < 103.8);
i3 = find(GridInfo.lat1 > 24);
i4 = find(GridInfo.lat1 < 1.5);
GridInfo.lon1(i1) = nan;
GridInfo.lon1(i2) = nan;
GridInfo.lon1(i3) = nan;
GridInfo.lon1(i4) = nan;
GridInfo.lat1(i1) = nan;
GridInfo.lat1(i2) = nan;
GridInfo.lat1(i3) = nan;
GridInfo.lat1(i4) = nan;
GridInfo.depth(i1) = nan;
GridInfo.depth(i2) = nan;
GridInfo.depth(i3) = nan;
GridInfo.depth(i4) = nan;
% for ndepth = [1:80:1000]
%     [C,h] = contour(GridInfo.lon1,GridInfo.lat1,-GridInfo.depth,[-ndepth 100]/10);
%     set(h, 'LineWidth',1.2)
%     % set(h, 'LineColor',[0.25 0 0.5])
%     set(h, 'LineColor','k')
%     % set(h, 'LineStyle','-.')
%     hh = get(h,'Children');
%     hold on
%     for i=1:numel(hh)
%         zdata = ones(size( get(hh(i),'XData') ));
%         set(hh(i), 'ZData',-ndepth*zdata/10+2)
%     end
% end
i1 = find(Coastal.lon<112);
Coastal.lon(i1) = [];
Coastal.lat(i1) = [];
i1 = find(Coastal.lat<21.5);
Coastal.lon(i1) = [];
Coastal.lat(i1) = [];
for nlevel = [1,5,9]%[1,7,13,19,24,30]
    FillLandMap_3D(Coastal.lon,Coastal.lat,1,nlevel);   % 画岸线及其陆地填充
    for i = [3,181]%1:size(GridInfo.lon2,1)
        h = plot(GridInfo.lon2(i,:),GridInfo.lat2(i,:),'Color',[0.75 0.75 0.75],'LineWidth',0.125);
        xdata = get(h,'XData');
        set(h,'ZData',-nlevel*ones(size(xdata)));
    end
    for i = [3,184]%1:size(GridInfo.lon2,2)
        h = plot(GridInfo.lon2(:,i),GridInfo.lat2(:,i),'Color',[0.75 0.75 0.75],'LineWidth',0.125);
        xdata = get(h,'XData');
        set(h,'ZData',-nlevel*ones(size(xdata)));
    end
end
caxis([1.5 2.3])
%legend('Layer 1 to 6','Layer 7 to 12','Layer 13 to 16')

subplot(1,2,2)
axis([112.5 115.6 20.5 23 -10 1])
view(25,30)
set(gca,'DataAspectRatio',[1.2 1 12/2.5])
hold on
load('G:\LB\碳循环\RCA_SFM_EOF\new\mycolor4.mat')
colormap(mycolor4);

for nobs = 1:size(Obs,1)
    add1 = floor(Station_Level(nobs)/6)*4;
    %     scatter3(GridInfo.lon(Station_Lon(nobs),Station_Lat(nobs)), ...
    %         GridInfo.lat(Station_Lon(nobs),Station_Lat(nobs)), ...
    %         -Station_Level(nobs)-add1, ...
    %         10,Observing(nobs),'filled');
    if Station_Lon(nobs) > 0
        switch floor(Station_Level(nobs)/6)
            case 0
                scatter3(GridInfo.lon(Station_Lon(nobs),Station_Lat(nobs)), ...
                    GridInfo.lat(Station_Lon(nobs),Station_Lat(nobs)), ...
                    -add1-1, ...
                    50,Simulating(nobs),'filled','^','MarkerEdgeColor','k');
                hold on
            case 1
                scatter3(GridInfo.lon(Station_Lon(nobs),Station_Lat(nobs)), ...
                    GridInfo.lat(Station_Lon(nobs),Station_Lat(nobs)), ...
                    -add1-1, ...
                    40,Simulating(nobs),'filled','MarkerEdgeColor','k');
                hold on
            case 2
                scatter3(GridInfo.lon(Station_Lon(nobs),Station_Lat(nobs)), ...
                    GridInfo.lat(Station_Lon(nobs),Station_Lat(nobs)), ...
                    -add1-1, ...
                    50,Simulating(nobs),'filled','v','MarkerEdgeColor','k');
                hold on
        end
    end
end
i1 = find(GridInfo.lon1 > 115);
i2 = find(GridInfo.lon1 < 103.8);
i3 = find(GridInfo.lat1 > 24);
i4 = find(GridInfo.lat1 < 1.5);
GridInfo.lon1(i1) = nan;
GridInfo.lon1(i2) = nan;
GridInfo.lon1(i3) = nan;
GridInfo.lon1(i4) = nan;
GridInfo.lat1(i1) = nan;
GridInfo.lat1(i2) = nan;
GridInfo.lat1(i3) = nan;
GridInfo.lat1(i4) = nan;
GridInfo.depth(i1) = nan;
GridInfo.depth(i2) = nan;
GridInfo.depth(i3) = nan;
GridInfo.depth(i4) = nan;
% for ndepth = [1:80:1000]
%     [C,h] = contour(GridInfo.lon1,GridInfo.lat1,-GridInfo.depth,[-ndepth 100]/10);
%     set(h, 'LineWidth',1.2)
%     % set(h, 'LineColor',[0.25 0 0.5])
%     set(h, 'LineColor','k')
%     % set(h, 'LineStyle','-.')
%     hh = get(h,'Children');
%     hold on
%     for i=1:numel(hh)
%         zdata = ones(size( get(hh(i),'XData') ));
%         set(hh(i), 'ZData',-ndepth*zdata/10+2)
%     end
% end
i1 = find(Coastal.lon<112);
Coastal.lon(i1) = [];
Coastal.lat(i1) = [];
i1 = find(Coastal.lat<21.5);
Coastal.lon(i1) = [];
Coastal.lat(i1) = [];
for nlevel = [1,5,9]%[1,7,13,19,24,30]
    FillLandMap_3D(Coastal.lon,Coastal.lat,1,nlevel);   % 画岸线及其陆地填充
    for i = [3,181]%1:size(GridInfo.lon2,1)
        h = plot(GridInfo.lon2(i,:),GridInfo.lat2(i,:),'Color',[0.75 0.75 0.75],'LineWidth',0.125);
        xdata = get(h,'XData');
        set(h,'ZData',-nlevel*ones(size(xdata)));
    end
    for i = [3,184]%1:size(GridInfo.lon2,2)
        h = plot(GridInfo.lon2(:,i),GridInfo.lat2(:,i),'Color',[0.75 0.75 0.75],'LineWidth',0.125);
        xdata = get(h,'XData');
        set(h,'ZData',-nlevel*ones(size(xdata)));
    end
end
caxis([1.5 2.3])
%legend('Layer 1 to 6','Layer 7 to 12','Layer 13 to 16')
