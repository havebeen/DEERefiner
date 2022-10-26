function [targetDistanceDistributionFinalForJSD, targetDistanceDistributionFinal] = targetDistanceDistributionFinalForJSDGenerator(distanceDistributionFullPath)
    for distanceDistributionNumber = 1:length(distanceDistributionFullPath)
        targetDistanceDistributionFinal(:, :, distanceDistributionNumber) = distanceDistributionFinalLoader(distanceDistributionFullPath(distanceDistributionNumber));
    end
    targetDistanceDistributionFinalForJSD = reshape(targetDistanceDistributionFinal(:, 2, :), [], 1);
end