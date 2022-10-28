function [distanceDistributionFinalForJSD, FINALPDBENSEMBLE] = distanceDistributionFinalForJSDGenerator(FINALPDBENSEMBLEFileName, residue1List, residue2List)
    FINALPDBENSEMBLE = pdbModelsLoader(FINALPDBENSEMBLEFileName);
    for FINALPDBNumber = 1:length(FINALPDBENSEMBLE(1, 1, :))
        for labelingIndex = 1:length(residue1List)
            [distanceDistributionFinal(:, labelingIndex, FINALPDBNumber), ~] = DEERefineMTSSLLabeling(FINALPDBENSEMBLE(:, :, FINALPDBNumber), residue1List(labelingIndex), residue2List(labelingIndex));
            distanceDistributionFinal(:, labelingIndex, FINALPDBNumber) = distanceDistributionFinal(:, labelingIndex, FINALPDBNumber)/sum(distanceDistributionFinal(:, labelingIndex, FINALPDBNumber));
        end
    end
    for FINALPDBNumber = 1:length(FINALPDBENSEMBLE(1, 1, :))
        distanceDistributionFinalForJSD(:, FINALPDBNumber) = reshape(distanceDistributionFinal(:, :, FINALPDBNumber), [], 1);
    end
end