function [movingCriterion, RMSE] = monteCarloMetropolisCriterionGenerator(oldDataSet, newDataSet, targetDataSet, temperature)
    oldDataSetNumbers = numel(oldDataSet);
    newDataSetNumbers = numel(newDataSet);
    targetDataSetNumbers = numel(targetDataSet);
    reshapedOldDataSet = reshape(oldDataSet, 1, oldDataSetNumbers);
    reshapedNewDataSet = reshape(newDataSet, 1, newDataSetNumbers);
    reshapedTargetDataSet = reshape(targetDataSet, 1, targetDataSetNumbers);
    differenceBetweenTargetAndOld = reshapedOldDataSet-reshapedTargetDataSet;
    differenceBetweenTargetAndNew = reshapedNewDataSet-reshapedTargetDataSet;
    oldRMSE = (sum(differenceBetweenTargetAndOld.^2)/oldDataSetNumbers)^0.5;
    newRMSE = (sum(differenceBetweenTargetAndNew.^2)/newDataSetNumbers)^0.5;
    RMSEDifference = oldRMSE - newRMSE;
    boltzmannFactor = exp(RMSEDifference/temperature);
    if RMSEDifference > 0
        movingCriterion = 1;
        RMSE = newRMSE;
    elseif RMSEDifference <= 0
        randomNumber = rand(1);
        if randomNumber > boltzmannFactor
            movingCriterion = 0;
            RMSE = oldRMSE;
        elseif randomNumber <= boltzmannFactor
            movingCriterion = 1;
            RMSE = newRMSE;
        end
    end
end