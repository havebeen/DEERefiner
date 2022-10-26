function candidateGenerator(candidateNumber, newCandidateFormatedStructurePhaseNine, timeStamp)
    finalPDBFileName = strcat("finalPDB_",  timeStamp, "_", num2str(candidateNumber), ".pdb");
    pdbSaver(finalPDBFileName, newCandidateFormatedStructurePhaseNine)
end