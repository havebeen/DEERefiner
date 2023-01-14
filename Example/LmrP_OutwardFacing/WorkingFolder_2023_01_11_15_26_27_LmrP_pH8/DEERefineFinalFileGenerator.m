function DEERefineFinalFileGenerator(residue1List, residue2List, FINALPDBENSEMBLE, minimalJSDIndex, targetDistanceDistributionFinal, distanceDistributionFullPath, timeStamp, R1_20210523)
    for labelingIndex = 1:length(residue1List)
        [DistanceDistributionFinal(:, 2, labelingIndex), DistanceDistributionFinal(:, 1, labelingIndex)] = DEERefineMTSSLLabeling(FINALPDBENSEMBLE(:, :, minimalJSDIndex), residue1List(labelingIndex), residue2List(labelingIndex), R1_20210523);
    end

    normalizedTargetDistanceDistributionFinal = targetDistanceDistributionFinal;
    for distanceDistributionNumber = 1:length(distanceDistributionFullPath)
        normalizedTargetDistanceDistributionFinal(:, 2, distanceDistributionNumber) = targetDistanceDistributionFinal(:, 2, distanceDistributionNumber)/max(targetDistanceDistributionFinal(:, 2, distanceDistributionNumber));
    end

    DEERefineFinalDistanceDistributionFileName = strcat("FINALPr_",  timeStamp, ".mat");
    DEERefineFinalDistanceDistribution = struct('targetDistanceDistributionFinal', normalizedTargetDistanceDistributionFinal,...
                                                'DistanceDistributionFinal', DistanceDistributionFinal, ...
                                                'residueIndex1', residue1List, ...
                                                'residueIndex2', residue2List);
    save(DEERefineFinalDistanceDistributionFileName, '-struct', 'DEERefineFinalDistanceDistribution');
    
    DEERefinerFigure = figure;
    figureHeight = 5;
    figureWiedth = ceil(distanceDistributionNumber/5);
    for distanceDistributionIndex = 1:distanceDistributionNumber
        subplot(figureHeight, figureWiedth, distanceDistributionIndex)
        hold on
        plot(normalizedTargetDistanceDistributionFinal(:, 1, distanceDistributionIndex)/10, normalizedTargetDistanceDistributionFinal(:, 2, distanceDistributionIndex), 'k', 'linewidth', 2.5);
        plot(DistanceDistributionFinal(:, 1, distanceDistributionIndex)/10, DistanceDistributionFinal(:, 2, distanceDistributionIndex), 'r', 'linewidth', 1.8);
        xlim([1 9])
        title(strcat(num2str(residue1List(distanceDistributionIndex)), '/', num2str(residue2List(distanceDistributionIndex)), 'R1'))
        xlabel('Distance, nm')
        box on
        set(gca,'FontSize',16)
        yticks([])
        ylabel('P(r)')
    end
    saveas(DEERefinerFigure, strcat("FINALPr_",  timeStamp, ".fig"))
    close(DEERefinerFigure)
end