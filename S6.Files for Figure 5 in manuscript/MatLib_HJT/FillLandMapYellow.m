function FillLandMapYellow(lon,lat)
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
       h=fill(lon(ns:IndNaN(i)-1),lat(ns:IndNaN(i)-1),[0.93 0.93 0.93]); 
       set(h,'linewidth',0.5); 
   else       
       plot(lon(ns:IndNaN(i)-1),lat(ns:IndNaN(i)-1),[0.93 0.93 0.93]);     
   end;   
   ns=IndNaN(i)+1;   
   if i== length(IndNaN)-1   
       h=fill(lon(ns:IndNaN(i+1)-1),lat(ns:IndNaN(i+1)-1),[0.93 0.93 0.93]);     
       set(h,'linewidth',0.5);   
   end     
   hold on; 
end; 
clear i
