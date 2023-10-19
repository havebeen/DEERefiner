function candidateGenerator(candidateNumber, newCandidateFormatedStructurePhaseNine, timeStamp, simulationIndex)
    finalPDBFileName = strcat("finalPDB_",  timeStamp, "_", num2str(simulationIndex), '_', num2str(candidateNumber), ".pdb");
    pdbSaver(finalPDBFileName, newCandidateFormatedStructurePhaseNine)
end