function Stats = calc_statistic(OBSTime,OBSVal,ModelTime,ModelVal)
% function Stats = calc_statistic(OBSTime,OBSVal,ModelTime,ModelVal)
% calculate statistic metrics
% Created by Jiatang, SYSU

bad = find(isnan(OBSVal));
OBSVal(bad) = [];
OBSTime(bad) = [];

NewOBSVal = [];
NewModelVal = [];
for ii = 1:numel(OBSVal)
    index = find(ModelTime == OBSTime(ii));
    if ~isempty(index)
        NewOBSVal = [NewOBSVal OBSVal(ii)];
        NewModelVal = [NewModelVal ModelVal(index)];
    end
end

corr0 = corrcoef(NewModelVal,NewOBSVal);
varD = NewOBSVal - mean(NewOBSVal);
RM = NewModelVal-NewOBSVal; % Model-data

Stats.RM = RM;
% Ranges of the observations and model results
Stats.RangeOBSVal = range(NewOBSVal);
Stats.RangeModelVal = range(NewModelVal);
Stats.RangeRatio = Stats.RangeModelVal/Stats.RangeOBSVal;  % Ratio of range (model/data)
% Standand deviation 
Stats.StdOBSVal = std(NewOBSVal);
Stats.StdModelVal = std(NewModelVal);
Stats.StdRatio = Stats.StdModelVal/Stats.StdOBSVal;  % Ratio of standard deviation
% Mean
Stats.MeanOBSVal = mean(NewOBSVal);
Stats.MeanModelVal = mean(NewModelVal);
% Bias (model-data)
Stats.Bias = Stats.MeanModelVal-Stats.MeanOBSVal;
Stats.Bias_RT_Mean = Stats.Bias/abs(Stats.MeanOBSVal+eps)*100;  % Bias relative to the observed mean
Stats.Bias_RT_Range = Stats.Bias/Stats.RangeOBSVal*100;  % Bias relative to the range of the observations
% RMS error
Stats.RMSE = sqrt(sum(RM.*RM)./length(RM));
Stats.RMSE_RT_Mean = Stats.RMSE/abs(Stats.MeanOBSVal+eps)*100;  % RMSE relative to the observed mean
Stats.RMSE_RT_Range = Stats.RMSE/Stats.RangeOBSVal*100;  % RMSE relative to the range of the observations
% Correlation
Stats.Corr = corr0(1,2);
% Model Efficiency
Stats.ME = 1-sum(RM.*RM)./sum(varD.*varD);

return
end
