function DEERefineFinalFileCopier(timeStamp)
    DEERefineFinalDistanceDistributionFileName = strcat("FINALPr_",  timeStamp, ".mat");
    FINALPDBENSEMBLEFileName = strcat("FINALPDBENSEMBLE_",  timeStamp, ".pdb");
    FINALPDBFileName = strcat("FINALPDB_",  timeStamp, ".pdb");
    copyfile(DEERefineFinalDistanceDistributionFileName, '..')
    copyfile(FINALPDBENSEMBLEFileName, '..')
    copyfile(FINALPDBFileName, '..')
end