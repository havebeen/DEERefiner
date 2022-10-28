function DEERefineFinalFileGenerator(residue1List, residue2List, FINALPDBENSEMBLE, minimalJSDIndex, targetDistanceDistributionFinal, distanceDistributionFullPath, timeStamp)
    for labelingIndex = 1:length(residue1List)
        [DistanceDistributionFinal(:, 2, labelingIndex), DistanceDistributionFinal(:, 1, labelingIndex)] = DEERefineMTSSLLabeling(FINALPDBENSEMBLE(:, :, minimalJSDIndex), residue1List(labelingIndex), residue2List(labelingIndex));
    end

    normalizedTargetDistanceDistributionFinal = targetDistanceDistributionFinal;
    for distanceDistributionNumber = 1:length(distanceDistributionFullPath)
        normalizedTargetDistanceDistributionFinal(:, 2, distanceDistributionNumber) = targetDistanceDistributionFinal(:, 2, distanceDistributionNumber)/max(targetDistanceDistributionFinal(:, 2, distanceDistributionNumber));
    end

    DEERefineFinalDistanceDistributionFileName = strcat("FINALPr_",  timeStamp, ".mat");
    DEERefineFinalDistanceDistribution = struct('targetDistanceDistributionFinal', normalizedTargetDistanceDistributionFinal,...
                                                'DistanceDistributionFinal', DistanceDistributionFinal);
    save(DEERefineFinalDistanceDistributionFileName, '-struct', 'DEERefineFinalDistanceDistribution');
end