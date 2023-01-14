function [distanceDistributionFinalForRMSE, FINALPDBENSEMBLE] = distanceDistributionFinalForJSDGenerator(FINALPDBENSEMBLEFileName, residue1List, residue2List, R1_20210523)
    FINALPDBENSEMBLE = pdbModelsLoader(FINALPDBENSEMBLEFileName);
    for FINALPDBNumber = 1:length(FINALPDBENSEMBLE(1, 1, :))
        for labelingIndex = 1:length(residue1List)
            [distanceDistributionFinal(:, labelingIndex, FINALPDBNumber), ~] = DEERefineMTSSLLabeling(FINALPDBENSEMBLE(:, :, FINALPDBNumber), residue1List(labelingIndex), residue2List(labelingIndex), R1_20210523);
        end
        distanceDistributionFinalForRMSE(:, FINALPDBNumber) = reshape(distanceDistributionFinal(:, :, FINALPDBNumber), [], 1);
    end
end