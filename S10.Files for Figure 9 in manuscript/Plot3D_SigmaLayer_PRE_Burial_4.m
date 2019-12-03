% 3D plot

clear all
close all
clc
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
GridInfo.XAZ = nan(GridInfo.nx,GridInfo.ny);
GridInfo.x = nan(GridInfo.nx,GridInfo.ny);
GridInfo.y = nan(GridInfo.nx,GridInfo.ny);
for irow = 1:length(griddata)
    jj = griddata(irow,1);  % 列
    ii = griddata(irow,2);  % 行
    GridInfo.lon(jj,ii) = griddata(irow,8);
    GridInfo.lat(jj,ii) = griddata(irow,7);
    GridInfo.depth(jj,ii) = griddata(irow,5);
    GridInfo.angle(jj,ii) = griddata(irow,6)*pi/180;
    GridInfo.x(jj,ii) = griddata(irow,3);
    GridInfo.y(jj,ii) = griddata(irow,4);
    GridInfo.XAZ(jj,ii) = griddata(irow,3)*griddata(irow,4);
    if GridInfo.depth(jj,ii) < 0 || GridInfo.depth(jj,ii) > 9999
        fprintf('%4d%4d%8.1f \n',jj,ii,GridInfo.depth(jj,ii));
        error('Too small or too large values found in depth!');
    else
        GridInfo.mask(jj,ii) = 1;
    end
end
clearvars griddata;
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
GridInfo.depth(isnan(GridInfo.depth)) = -0.0;
[NewLon,NewLat] = meshgrid(112:0.004:116,20:0.004:23);
GridInfo.lon(isnan(GridInfo.lon)) = -0.0;
GridInfo.lat(isnan(GridInfo.lat)) = -0.0;
WetGrid.Index = find(~isnan(GridInfo.lon));
WetGrid.Lon = GridInfo.lon(WetGrid.Index);
WetGrid.Lat = GridInfo.lat(WetGrid.Index);
load('F:\LB\碳循环\RCA_SFM_EOF\AreaAndVol.mat');
area = AreaAndVol.area;
vol = AreaAndVol.vol;
h = vol./repmat(area,[1,1,size(vol,3)]);
    figure(1)
    
    set(gca,'DataAspectRatio',[4/5 1 100])
    
    axis([113 115 21.5 23 -220 1])
    hold on
    load('F:\LB\碳循环\RCA_SFM_EOF\new\mycolor4.mat')
    colormap(mycolor4);

nmon = {'Jan','Feb','Mar','Apr','May','Jun', ...
    'Jul','Aug','Sep','Oct','Nov','Nov'};

for imonth = 1:2:12
    clear Data
    AO1 = [];
    AO2 = [];
    AO3 = [];
    AO4 = [];
    load(['F:\LB\碳循环\RCA_SFM_EOF\3DPlot\monthlyVar\',nmon{imonth},'_JPOCN.mat']);
    AO1(:,:) = mean(JPOCN,4).*area;
    load(['F:\LB\碳循环\RCA_SFM_EOF\3DPlot\monthlyVar\',nmon{imonth},'_JDIC.mat']);
    AO2(:,:) = mean(JDIC,4);
    Data = AO1/12 -1000*(AO2);
    
    Data(3,:) = nan;
    Data(:,3) = nan;
    Data(181,:) = nan;
    Data = repmat(Data,[1 1 16]);
    % Data(find(Data>10000)) = nan;
    % Data(find(Data <= 0)) = nan;
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
                    Xnew(ii,jj,1:abs(round(max((Z1(1:end-1)))*10))+1) = Xnew1(i1(end));
                end
            end
            
            if GridInfo.depth(ii,jj) > 0
                Z2 = [GridInfo.depth(ii,jj)];
                index1_Z2 = round(Z2(1)*10);
                Z3(ii,jj,[index1_Z2+1:1000+1]) = [-1*index1_Z2/10:-0.1:-100];
            else
                Z3(ii,jj,1:end) = [0:-0.1:-100];
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
    V = Xnew;

    for ndepth = imonth%[1,50,100,150,200]
        NewLon1 = NewLon;
        NewLat1 = NewLon;
        if 1
            OriginalModelData = V(:,:,1);
            GoodData = OriginalModelData(WetGrid.Index);
            F = TriScatteredInterp(WetGrid.Lon,WetGrid.Lat,GoodData);
            InterpData = F(NewLon,NewLat);
            i1 = find(NewLon > 114.1);
            i2 = find(NewLon < 112.9);
            i3 = find(NewLat > 22.8);
            i4 = find(NewLat < 21.8);
            i1 = find(NewLon > 115);
            i2 = find(NewLon < 113);
            i3 = find(NewLat > 25);
            i4 = find(NewLat < 21.8);
            InterpData(i1) = nan;
            InterpData(i2) = nan;
            InterpData(i3) = nan;
            InterpData(i4) = nan;
            NewLon1(i1) = nan;
            NewLon1(i2) = nan;
            NewLon1(i3) = nan;
            NewLon1(i4) = nan;
            NewLat1(i1) = nan;
            NewLat1(i2) = nan;
            NewLat1(i3) = nan;
            NewLat1(i4) = nan;
            PropertyOpt.Colors{1} = 7.1:0.0125:8.3;
            PropertyOpt.Colors{1} = 1.3:0.0125:2.1;
            [h] = pcolor(NewLon,NewLat-0.25*(ndepth-1)/2,InterpData); hold on
            shading flat
            shading interp
            zdata = ones(size(get(h,'XData')));
            set(h,'ZData',-ndepth*zdata*25)
            
            [C,h] = contour(NewLon,NewLat-0.25*(ndepth-1)/2,InterpData,[0 0]); hold on
            hh = get(h,'Children');
            for ih = 1:numel(hh)
                zdata = ones(size(get(hh(ih),'XData')));
                set(hh(ih),'ZData',-ndepth*zdata*25)
            end
            set(h,'Color','w')
            set(h,'LineWidth',1)
        end
    end
end
i1 = find(Coastal.lat < 20);
Coastal.lat(i1) = [];
Coastal.lon(i1) = [];
for ndepth = 1:2:12
    FillLandMap_3D(Coastal.lon,Coastal.lat-0.25*(ndepth-1)/2,1,ndepth*25-1);
    hold on
end
colorbar
view(30,30)
caxis([-100 100])
 axis([113 115 19.5 23 -220 1])
% colorbar