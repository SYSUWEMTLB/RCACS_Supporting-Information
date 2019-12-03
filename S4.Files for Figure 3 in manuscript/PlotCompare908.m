% How to create a taylor diagram with overlaid markers
%
% An eighth example of how to create a Taylor diagram given one set
% of reference observations and multiple model predictions for the
% quantity.
%
% This example is a variation on the seventh example (taylor7) where now a
% fourth data point having a negative correlation is overlaid on an
% existing Taylor diagram that already has 3 data points with positive
% correlations. It is chosen to have data points with positive correlations
% appear in red while data points with negative correlations are displayed
% in blue.
%
% All functions in the Skill Metrics Toolbox are designed to only work with
% one-dimensional arrays, e.g. time series of observations at a selected
% location. The one-dimensional data are read in as data structures via a
% mat file. The latter are stored in data structures in the format:
% ref.data, pred1.data, pred2.dat, and pred3.dat. The plot is written to a
% file in Portable Network Graphics (PNG) format.
%
% The reference data used in this example are cell concentrations of a
% phytoplankton collected from cruise surveys at selected locations and
% time. The model predictions are from three different simulations that
% have been space-time interpolated to the location and time of the sample
% collection. Details on the contents of the data structures (once loaded)
% can be obtained by simply entering the data structure variable name at
% the command prompt, e.g.
% >> ref
% ref =
%          data: [57x1 double]
%          date: {57x1 cell}
%         depth: [57x1 double]
%      latitude: [57x1 double]
%     longitude: [57x1 double]
%       station: [57x1 double]
%          time: {57x1 cell}
%         units: 'cell/L'
%          jday: [57x1 double]

% Author: Peter A. Rochford
%         Symplectic, LLC
%         www.thesymplectic.com
%         prochford@thesymplectic.com

% Close any previously open graphics windows
clear all
clc
close all;
clf
addpath('G:\LB\Ì¼Ñ­»·\RCAcarboncycle\PeterRochford-SkillMetricsToolbox-7ceead9');
% Set the figure properties (optional)
%set(gcf,'units','inches','position',[0,10.0,14.0,10.0]);
set(gcf,'units','inches','position',[1,-3,14.0,10.0]);
set(gcf, 'DefaultLineLineWidth', 1.5); % linewidth for plots
set(gcf,'DefaultAxesFontSize',18); % font size of axes text

% Read in data from a mat file
load('Compare_06Jul_1202.mat')
if 1
    Var = {' ecomtemp',' salt','talk','pH','pCO2','dic'};
    name = {' Temp',' Sal',' TALK',' pH',' pCO_2',' DIC'};
    color = {[255/255 0 0],[255/255 0 0], ...
        [0.11765 0.56471 1],[0.11765 0.56471 1],[0.11765 0.56471 1],[0.11765 0.56471 1]
        };
    kind={'o','o','s','s','s','s','s'};
end
for ii = 1:numel(Var)
    eval(['pred.data = Compare.',Var{ii},'(:,2);']);
    eval(['ref.data = Compare.',Var{ii},'(:,1);']);
    eval(['taylor_stats',num2str(ii),' = taylor_statistics(pred,ref,''data'');'])
    clear pred ref
    
    eval(['sdev = [1 ; taylor_stats',num2str(ii),'.sdev(2)/taylor_stats',num2str(ii),'.sdev(1)];']);
    eval(['crmsd = [0 ; taylor_stats',num2str(ii),'.crmsd(2)];']);
    eval(['ccoef = [1 ; taylor_stats',num2str(ii),'.ccoef(2)];']);
    sdev1 = [1;1];
    crmsd1 = [0;0];
    ccoef1 = [1;1];
    label = Var{ii};
    if ii == 1
        [hp, ht, axl] = taylor_diagram(sdev,crmsd,ccoef, ...
            'markerLabelColor', color{ii}, 'markerColor',color{ii}, ...
            'tickRMS',0.0:0.3:1.5,'markerLabel',name{ii}, ...
            'colRMS',[105/255 205/255 111/255], 'styleRMS', ':', 'widthRMS', 2.0, 'titleRMS', 'on', ...
            'tickSTD',0.0:0.5:2, 'limSTD',2, ...
            'colSTD','b', 'styleSTD', '-.', 'widthSTD', 2, 'titleSTD', 'on', ...
            'colCOR','k', 'styleCOR', '--', 'widthCOR', 1.0, 'titleCOR', 'on', ...
            'markerSize',12, 'markerObs','p','colOBS',[255/255 0/255 0/255],...
            'markerKind',kind{ii},'MarkerFaceColor',color{ii});
    else
        taylor_diagram(sdev,crmsd,ccoef, ...
            'overlay','on','markerLabel',name{ii}, ...
            'markerLabelColor', color{ii}, ...
            'markerColor',color{ii},'markerKind',kind{ii},'MarkerFaceColor',color{ii});
%         taylor_diagram(sdev,crmsd,ccoef, ...
%             'overlay','on', ...
%             'markerLabelColor', color{ii}, ...
%             'markerColor',color{ii},'markerKind',kind{ii},'MarkerFaceColor',color{ii});
    end
    clear sdev crmsd ccoef
end
writepng(gcf,'Compare908.png');