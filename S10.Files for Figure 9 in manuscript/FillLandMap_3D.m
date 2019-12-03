function FillLandMap_3D(lon,lat,index,index1)
%FillLandMap fill land points
%
% Example:
%     FillLandMap(lon,lat)
%
% Modified by Jiatang Hu on 2013/05/17

ns=1;
IndNaN=find(isnan(lon));
for i=1:290 %length(IndNaN)-1
    if IndNaN(i+1)-IndNaN(i)>2
        %        h=fill(lon(ns:IndNaN(i)-1),lat(ns:IndNaN(i)-1),'w');
        h=fill3(lon(ns:IndNaN(i)-1),lat(ns:IndNaN(i)-1),-1*index1*ones(size(lon(ns:IndNaN(i)-1))),[1 1 1]*index);
        set(h,'linewidth',1);

%             set(h, 'ZData',index1)

    else
        plot(lon(ns:IndNaN(i)-1),lat(ns:IndNaN(i)-1),'k');
   %     set(gca,'ZData',index1);
    end;
    ns=IndNaN(i)+1;
    if i== length(IndNaN)-1
        h=fill3(lon(ns:IndNaN(i+1)-1),lat(ns:IndNaN(i+1)-1),-1*index1*ones(size(lon(ns:IndNaN(i)-1))),'k');
        set(h,'linewidth',1);
%             set(h, 'ZData',index1)
    end
    hold on;
end;
clear i
