function LOC = FindModelGridPoints(Lon,Lat,gridfile)

% Load grid data:
griddata = load(gridfile);

AllIX = griddata(:,1);
AllIY = griddata(:,2);
AllLat = griddata(:,7);
AllLon = griddata(:,8);
    
rad = pi/180;

DiffLon = abs(AllLon-Lon);
DiffLat = abs(AllLat-Lat);
% Distant = sqrt((DiffLon.*cos(AllLat.*rad)).^2+DiffLat.^2);
Distant = sqrt(DiffLon.^2+DiffLat.^2);

[DMin,Index] = min(Distant);
LOC.ix = AllIX(Index);
LOC.iy = AllIY(Index);

LOC.Bdep = griddata(Index,5);

return
end

