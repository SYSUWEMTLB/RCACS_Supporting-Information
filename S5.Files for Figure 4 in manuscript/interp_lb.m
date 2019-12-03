
function [ib] = interp_lb(b)
s=size(b);
ind=find(~isnan(b));
[i j]=ind2sub(s,ind);
v=b(ind);
[ii jj]=ndgrid(1:s(1),1:s(2));
ib=griddata(i,j,v,ii,jj);
end
