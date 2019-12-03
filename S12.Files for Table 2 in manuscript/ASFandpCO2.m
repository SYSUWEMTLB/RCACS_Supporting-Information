clear all
close all
clc

load('matfile\Feb_pCO2.mat');
load('matfile\Feb_ASF.mat');

nregion = 4; % 分区数量
for nr = 1:nregion
    switch nr
        case 1
            c = [30,30,65,65]; % 划分点列号
            r = [40,60,60,40]; % 划分点行号
            i1 = [1,2,3,4,1]; % 区域1划分点序号 需要回到出发点形成闭环
        case 2
            c = [90,90,183,183]; % 划分点列号
            r = [70,129,129,70]; % 划分点行号
            i2 = [1,2,3,4,1]; % 区域1划分点序号 需要回到出发点形成闭环
        case 3
            c = [55,55,89,89]; % 划分点列号
            r = [71,130,130,71]; % 划分点行号
            i3 = [1,2,3,4,1]; % 区域1划分点序号 需要回到出发点形成闭环
        case 4
            c = [30,30,54,54]; % 划分点列号
            r = [61,130,130,61]; % 划分点行号
            i4 = [1,2,3,4,1]; % 区域1划分点序号 需要回到出发点形成闭环
    end
    clear BW cnr rnr
    pCO2 = mean(pCO2,4);
    pCO2(find(pCO2 == 0)) = nan;
    pCO2 = pCO2(:,:,1);
    eval(['cnr = c(i',num2str(nr),');'])
    eval(['rnr = r(i',num2str(nr),');'])
    [BW,xi,yi] = roipoly(pCO2,cnr,rnr);
    eval(['Region_pCO2.region',num2str(nr),' = BW.*pCO2;']);
    
    clear BW cnr rnr
    ASF = mean(ASF,4);
    ASF(find(ASF == 0)) = nan;
    ASF = ASF(:,:,1);
    eval(['cnr = c(i',num2str(nr),');'])
    eval(['rnr = r(i',num2str(nr),');'])
    [BW,xi,yi] = roipoly(ASF,cnr,rnr);
    eval(['Region_ASF.region',num2str(nr),' = BW.*ASF;']);
end

min(Region_pCO2.region1(find(abs(Region_pCO2.region1)>0)))
max(Region_pCO2.region1(find(abs(Region_pCO2.region1)>0)))
mean(Region_ASF.region1(find(abs(Region_ASF.region1)>0)))
std(Region_ASF.region1(find(abs(Region_ASF.region1)>0)))

min(Region_pCO2.region2(find(abs(Region_pCO2.region2)>0)))
max(Region_pCO2.region2(find(abs(Region_pCO2.region2)>0)))
mean(Region_ASF.region2(find(abs(Region_ASF.region2)>0)))
std(Region_ASF.region2(find(abs(Region_ASF.region2)>0)))

min(Region_pCO2.region3(find(abs(Region_pCO2.region3)>0)))
max(Region_pCO2.region3(find(abs(Region_pCO2.region3)>0)))
mean(Region_ASF.region3(find(abs(Region_ASF.region3)>0)))
std(Region_ASF.region3(find(abs(Region_ASF.region3)>0)))

min(Region_pCO2.region4(find(abs(Region_pCO2.region4)>0)))
max(Region_pCO2.region4(find(abs(Region_pCO2.region4)>0)))
mean(Region_ASF.region4(find(abs(Region_ASF.region4)>0)))
std(Region_ASF.region4(find(abs(Region_ASF.region4)>0)))

