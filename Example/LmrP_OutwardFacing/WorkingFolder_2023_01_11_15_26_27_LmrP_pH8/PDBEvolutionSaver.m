function PDBEvolutionSaver(outputTrajectoryPDBFileName, formatedStructure, PDBEvolution)
    if PDBEvolution == 0
        return
    elseif PDBEvolution == 1
        pdbSaver(outputTrajectoryPDBFileName, formatedStructure, 1, 1);
    end
end