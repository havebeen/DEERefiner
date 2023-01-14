function FINALPDBENSEMBLEFileGenerator(timeStamp)
    FINALPDBENSEMBLEFileName = strcat("FINALPDBENSEMBLE_",  timeStamp, ".pdb");
    finalPDBFileNameHeader = strcat("finalPDB_",  timeStamp);
    finalPDBFiles = dir(strcat('*', finalPDBFileNameHeader, '*'));
    finalPDBNumber = length(finalPDBFiles);
    for loopFinalPDBNumber = 1:finalPDBNumber
        FINALPDBFileNametmp = finalPDBFiles(loopFinalPDBNumber).name;
        FINALPDBtmp(:, :, loopFinalPDBNumber) = pdbLoader(FINALPDBFileNametmp);
        pdbSaver(FINALPDBENSEMBLEFileName, FINALPDBtmp(:, :, loopFinalPDBNumber), 1, loopFinalPDBNumber);
    end
end