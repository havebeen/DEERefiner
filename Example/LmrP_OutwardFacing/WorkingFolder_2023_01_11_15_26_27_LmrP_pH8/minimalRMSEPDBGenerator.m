function minimalRMSEIndex = minimalRMSEPDBGenerator(FINALPDBENSEMBLE, FINALPDBMaximalDifference, minimalFINALPDBMaximalDifference, FINALPDBRMSE, minimalFINALPDBRMSE, timeStamp)
    FINALPDBFileName = strcat("FINALPDB_",  timeStamp, ".pdb");
    minimalRMSEIndex = find((FINALPDBMaximalDifference == minimalFINALPDBMaximalDifference &FINALPDBRMSE== minimalFINALPDBRMSE)==1);
    pdbSaver(FINALPDBFileName, FINALPDBENSEMBLE(:, :, minimalRMSEIndex(1)));
end