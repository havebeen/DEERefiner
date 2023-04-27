R1LibraryMatFileName = 'R1_library_20210523.mat';
load(R1LibraryMatFileName);
DEERefineSturctureFileName = "DEERefineSturcture.mat";
load(DEERefineSturctureFileName);
structureIndexEnd = structureIndexEachScript;

for structureIndex = 1:structureIndexEnd
    clearvars -except structureIndex simulationIndex structureIndexEnd app ...
                      runFileName RMSE maximalClashes monteCarloIteration initialStructureFullPath phiPsiAngle flexiblePhiPsiIndex targetNumberStructure...
                      R1LibraryMatFileName R1_20210523 PDBEvolution monteCarloTemperature
                  
    DEERefineSturctureFileName = "DEERefineSturcture.mat";
    load(DEERefineSturctureFileName);
    rng('shuffle')
    timeStamp = char(erase(runFileName(1),"runFile_"));
    timeStamp(end-1:end) = [];
    currentRunFileName = string(mfilename);
    simulationIndex = erase(currentRunFileName, strcat("runFile_", timeStamp, "_"));
    currentStatusFileName = strcat("currentStatus_", timeStamp, ".dat");
    advancedPrFileName = strcat(timeStamp, '_advanced_Pr.dat');
    firstRMSEPassedNumberFileName = strcat(timeStamp, '_first_RMSE_passed_number.dat');
    secondRMSEPassedNumberFileName = strcat(timeStamp, '_second_RMSE_passed_number.dat');
    thirdRMSEPassedNumberFileName = strcat(timeStamp, '_third_RMSE_passed_number.dat');
    fourthRMSEPassedNumberFileName = strcat(timeStamp, '_fourth_RMSE_passed_number.dat');
    fifthRMSEPassedNumberFileName = strcat(timeStamp, '_fifth_RMSE_passed_number.dat');
    
    numberOfStructure = targetNumberStructure;
    fifthRMSEPassedNumber = load(fifthRMSEPassedNumberFileName);
    if fifthRMSEPassedNumber >= numberOfStructure
        return
    end
%     maximalClashes = app.maximalClashes;
    monteCarloIterationNumber = monteCarloIteration;
%     RMSE = app.RMSE;
    [residue1List, residue2List, distanceDistributionList] = advancedPrFileLoader(advancedPrFileName);  % Distance distribution files loading
    
    initialFullFormatedPDB = pdbLoader(initialStructureFullPath);   % initial structure pdb
    
    initialMutatedFormatedPDB = brokenFormatedPDBFixer(initialFullFormatedPDB); % mutated structure pdb
    
    sideChainRotatorResidueIndex = unique(double(initialMutatedFormatedPDB(:, 7))); %% side chain rotator residue index
    
    initialFormatedBackbone = initialMutatedFormatedPDB(initialMutatedFormatedPDB(:, 3) == " N  "|...
                                                        initialMutatedFormatedPDB(:, 3) == " CA "|...
                                                        initialMutatedFormatedPDB(:, 3) == " C  "|...
                                                        initialMutatedFormatedPDB(:, 3) == " O  ", :);  % extracting initial formated backbone to obtain backbone coordinates for structural alignment and calculating minimalBackboneNonBondedDistance.
    
    initialBackboneGeometry = formatedBackbone2Geometry(initialFormatedBackbone);
    
    initialBackboneCoordinates = double(initialFormatedBackbone(:, 9:11));  % initial backbone coordinates
    
    minimalBackboneNonBondedDistance = formatedBackbone2MinimalBackboneNonBondedDistance(initialFormatedBackbone); % backbone constraints
    
    
    for MTSSLLabelingIndex = 1:length(residue1List)
        [initialSimulatedDistanceDistributionY(:, MTSSLLabelingIndex), ~] = DEERefineMTSSLLabeling(initialMutatedFormatedPDB, residue1List(MTSSLLabelingIndex), residue2List(MTSSLLabelingIndex), R1_20210523);
    end
    
    outputTrajectoryPDBFileName = strcat(timeStamp, '_Trajectory_', num2str(simulationIndex), '_', num2str(structureIndex), '.pdb');
    
    normalizedInitialSimulatedDistanceDistributionY = initialSimulatedDistanceDistributionY./sum(initialSimulatedDistanceDistributionY);
    
    [~, targetDEERMaximalIndexTmp] = max(distanceDistributionList);
    
    targetDEERMaximalIndex = reshape(targetDEERMaximalIndexTmp(1, 2, :), 1, []);
    
    [~, initialDEERMaximalIndex] = max(normalizedInitialSimulatedDistanceDistributionY);
    
    RMSEFirstPhase(1) = (mean((targetDEERMaximalIndex-initialDEERMaximalIndex).^2))^0.5;
    
    iterationIndexFirstPhase(1) = 1;
    
    outputTrajectoryDatFileName = strcat(timeStamp, '_Trajectory_', num2str(simulationIndex), '_', num2str(structureIndex), '.dat');
    
    phiPsiAngleVariationAtLoop = phiPsiAngle;
    
    allResidueIndex = unique(double(initialMutatedFormatedPDB(:, 7)));

    flexibleRegionIndex = flexiblePhiPsiIndex;
    
    monteCarloOldGeometry = initialBackboneGeometry;
    monteCarloOldDistanceDistribution = initialDEERMaximalIndex;
    monteCarloOldStructure = initialMutatedFormatedPDB;
    monteCarloOldContactedResidueIndex = formatedPDB2contactedResidueNumbers(monteCarloOldStructure, 2.4);
    
    clashedResidueNumber(1) = length(monteCarloOldContactedResidueIndex);
    
    for monteCarloSteps = 1:monteCarloIterationNumber
        monteCarloNewCandidateGeometry = monteCarloOldGeometry;
        monteCarloNewCandidateGeometry.dihedral.CN1CA1C1(flexibleRegionIndex) = monteCarloNewCandidateGeometry.dihedral.CN1CA1C1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoop;
        monteCarloNewCandidateGeometry.dihedral.NCACN1(flexibleRegionIndex) = monteCarloNewCandidateGeometry.dihedral.NCACN1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoop;
        currentCandidateBackboneCoordinates = geometry2backboneCoordinates(monteCarloNewCandidateGeometry);
        [~, currentCandidateAlignedBackboneCoordinates] = procrustes(initialBackboneCoordinates, currentCandidateBackboneCoordinates, 'scaling', false);
        currentCandidateFormatedBackbone = initialFormatedBackbone;
        currentCandidateFormatedBackbone(:, 9:11) = string(currentCandidateAlignedBackboneCoordinates);
%         currentCandidateMinimalBackboneNonBondedDistance = formatedBackbone2MinimalBackboneNonBondedDistance(currentCandidateFormatedBackbone);
%         if minimalBackboneNonBondedDistance*0.98 > currentCandidateMinimalBackboneNonBondedDistance 
%             currentStatusUpdator(currentStatusFileName, 1)
%             continue
%         end
        currentCandidateFormatedStructure = oldSideChainInstaller(currentCandidateFormatedBackbone, monteCarloOldStructure);
        newCandidateFormatedStructure = sideChainRotator(currentCandidateFormatedStructure, 10, sideChainRotatorResidueIndex);
        newCandidateConteactedResidueIndex = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructure, 2.4);
        
        for labelingIndex = 1:length(residue1List)
            [newCandidateSimulatedDistanceDistributionY(:, labelingIndex), ~] = DEERefineMTSSLLabeling(newCandidateFormatedStructure, residue1List(labelingIndex), residue2List(labelingIndex), R1_20210523);
        end
        [~, newCandidateDEERMaximalIndex] = max(newCandidateSimulatedDistanceDistributionY);
        monteCarloNewDistanceDistribution = newCandidateDEERMaximalIndex;
        monteCarloClashesTemperature = 3;
        [DEERIfMovingCriterion, RMSEFirstPhase(monteCarloSteps+1)] = monteCarloMetropolisCriterionGenerator(monteCarloOldDistanceDistribution, monteCarloNewDistanceDistribution, targetDEERMaximalIndex, monteCarloTemperature);
        [clashesIfMovingCriterion, clashedResidueNumber(monteCarloSteps+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndex), length(newCandidateConteactedResidueIndex), 0, monteCarloClashesTemperature);
        if DEERIfMovingCriterion == 1 && clashesIfMovingCriterion == 1
            currentStatusUpdator(currentStatusFileName, 2)
            PDBEvolutionSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructure, PDBEvolution)
%             RMSEDatSaver(outputTrajectoryDatFileName, monteCarloSteps, RMSEFirstPhase(monteCarloSteps+1))
            datSaver(outputTrajectoryDatFileName, monteCarloSteps, RMSEFirstPhase(monteCarloSteps+1), clashedResidueNumber(monteCarloSteps+1))
            monteCarloOldContactedResidueIndex = newCandidateConteactedResidueIndex;
            monteCarloOldStructure = newCandidateFormatedStructure;
            monteCarloOldDistanceDistribution = monteCarloNewDistanceDistribution;
            monteCarloOldGeometry = monteCarloNewCandidateGeometry;
            if RMSEFirstPhase(monteCarloSteps+1) <= RMSE

                increaseFirstRMSEPassedNumber(firstRMSEPassedNumberFileName)
                currentStatusUpdator(currentStatusFileName, 3)

                monteCarloOldGeometryPhaseTwo = monteCarloOldGeometry;
                monteCarloOldStructurePhaseTwo = monteCarloOldStructure;
    %             pdbSaver(phaseOnePassedFileName, newCandidateFormatedStructure);
                break
            end
        else
            monteCarloOldContactedResidueIndex = monteCarloOldContactedResidueIndex;
            monteCarloOldStructure = monteCarloOldStructure;
            monteCarloOldDistanceDistribution = monteCarloOldDistanceDistribution;
            monteCarloOldGeometry = monteCarloOldGeometry;
        end
%         phaseOnePassedFileName = strcat(timeStamp, "_phase1_passed", num2str(simulationIndex), '_', num2str(structureIndex), '.pdb');
    end
    if exist('monteCarloOldStructurePhaseTwo', 'var') ~= 1
        continue
    end
    clearvars -except simulationIndex structureIndex structureIndexEnd timeStamp...
        currentStatusFileName...
        firstRMSEPassedNumberFileName secondRMSEPassedNumberFileName thirdRMSEPassedNumberFileName fourthRMSEPassedNumberFileName fifthRMSEPassedNumberFileName...
        residue1List residue2List distanceDistributionList...
        initialBackboneCoordinates initialFormatedBackbone...
        phiPsiAngleVariationAtLoop flexibleRegionIndex...
        sideChainRotatorResidueIndex minimalBackboneNonBondedDistance...
        monteCarloOldGeometryPhaseTwo monteCarloOldStructurePhaseTwo...
        maximalClashes RMSE app numberOfStructure runFileName distanceDistributionFullPath...
                      R1LibraryMatFileName R1_20210523 PDBEvolution...
                      outputTrajectoryPDBFileName outputTrajectoryDatFileName monteCarloTemperature monteCarloIterationNumber
        
    
    phiPsiAngleVariationAtLoopPhaseTwo = 0.05;
    clashesCriterionPhaseTwo = 2.3;
    
    
    monteCarloOldContactedResidueIndexPhaseTwo = formatedPDB2contactedResidueNumbers(monteCarloOldStructurePhaseTwo, clashesCriterionPhaseTwo);
    clashedResidueNumberPhaseTwo(1) = length(monteCarloOldContactedResidueIndexPhaseTwo);
    for phaseTwoSteps = 1:monteCarloIterationNumber
        fifthRMSEPassedNumber = load(fifthRMSEPassedNumberFileName);
        if fifthRMSEPassedNumber >= numberOfStructure
            return
        end
        monteCarloNewCandidateGeometryPhaseTwo = monteCarloOldGeometryPhaseTwo;
        monteCarloNewCandidateGeometryPhaseTwo.dihedral.CN1CA1C1(flexibleRegionIndex) = monteCarloNewCandidateGeometryPhaseTwo.dihedral.CN1CA1C1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoopPhaseTwo;
        monteCarloNewCandidateGeometryPhaseTwo.dihedral.NCACN1(flexibleRegionIndex) = monteCarloNewCandidateGeometryPhaseTwo.dihedral.NCACN1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoopPhaseTwo;
        currentCandidateBackboneCoordinatesPhaseTwo = geometry2backboneCoordinates(monteCarloNewCandidateGeometryPhaseTwo);
        [~, currentCandidateAlignedBackboneCoordinatesPhaseTwo] = procrustes(initialBackboneCoordinates, currentCandidateBackboneCoordinatesPhaseTwo, 'scaling', false);
        currentCandidateFormatedBackbonePhaseTwo = initialFormatedBackbone;
        currentCandidateFormatedBackbonePhaseTwo(:, 9:11) = string(currentCandidateAlignedBackboneCoordinatesPhaseTwo);
%         currentCandidateMinimalBackboneNonBondedDistance = formatedBackbone2MinimalBackboneNonBondedDistance(currentCandidateFormatedBackbonePhaseTwo);
%         if minimalBackboneNonBondedDistance*0.98 > currentCandidateMinimalBackboneNonBondedDistance 
%             continue
%         end
        currentCandidateFormatedStructurePhaseTwo = oldSideChainInstaller(currentCandidateFormatedBackbonePhaseTwo, monteCarloOldStructurePhaseTwo);
        newCandidateFormatedStructurePhaseTwo = sideChainRotator(currentCandidateFormatedStructurePhaseTwo, 25, monteCarloOldContactedResidueIndexPhaseTwo);
        newCandidateConteactedResidueIndexPhaseTwo = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructurePhaseTwo, clashesCriterionPhaseTwo);
        [clashesIfMovingCriterionPhaseTwo, clashedResidueNumberPhaseTwo(phaseTwoSteps+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndexPhaseTwo), length(newCandidateConteactedResidueIndexPhaseTwo), 0, 0.6);
        if clashesIfMovingCriterionPhaseTwo==1
            currentStatusUpdator(currentStatusFileName, 4)
            monteCarloOldContactedResidueIndexPhaseTwo = newCandidateConteactedResidueIndexPhaseTwo;
            monteCarloOldStructurePhaseTwo = newCandidateFormatedStructurePhaseTwo;
            monteCarloOldGeometryPhaseTwo = monteCarloNewCandidateGeometryPhaseTwo;
            PDBEvolutionSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructurePhaseTwo, PDBEvolution)
            datSaver(outputTrajectoryDatFileName, phaseTwoSteps, 0, clashedResidueNumberPhaseTwo(phaseTwoSteps+1))
            if clashedResidueNumberPhaseTwo(phaseTwoSteps+1) <= maximalClashes
                currentStatusUpdator(currentStatusFileName, 5)
                monteCarloOldGeometryPhaseThree = monteCarloNewCandidateGeometryPhaseTwo;
                monteCarloOldStructurePhaseThree = newCandidateFormatedStructurePhaseTwo;
                break
            end 
        else
            monteCarloOldContactedResidueIndexPhaseTwo = monteCarloOldContactedResidueIndexPhaseTwo;
            monteCarloOldStructurePhaseTwo = monteCarloOldStructurePhaseTwo;
            monteCarloOldGeometryPhaseTwo = monteCarloOldGeometryPhaseTwo;
        end
    end
    if exist('monteCarloOldStructurePhaseThree', 'var') ~= 1
        continue
    end
    clearvars -except simulationIndex structureIndex structureIndexEnd timeStamp...
        currentStatusFileName...
        firstRMSEPassedNumberFileName secondRMSEPassedNumberFileName thirdRMSEPassedNumberFileName fourthRMSEPassedNumberFileName fifthRMSEPassedNumberFileName...
        residue1List residue2List distanceDistributionList...
        initialBackboneCoordinates initialFormatedBackbone...
        sideChainRotatorResidueIndex minimalBackboneNonBondedDistance...
        phiPsiAngleVariationAtLoop flexibleRegionIndex...
        monteCarloOldGeometryPhaseThree monteCarloOldStructurePhaseThree...
        maximalClashes RMSE app numberOfStructure runFileName distanceDistributionFullPath...
                      R1LibraryMatFileName R1_20210523 PDBEvolution...
                      outputTrajectoryPDBFileName outputTrajectoryDatFileName monteCarloTemperature monteCarloIterationNumber
    phiPsiAngleVariationAtLoopPhaseThree = 0.05;
    
    for MTSSLLabelingIndex = 1:length(residue1List)
        [initialSimulatedDistanceDistributionYPhaseThree(:, MTSSLLabelingIndex), ~] = DEERefineMTSSLLabeling(monteCarloOldStructurePhaseThree, residue1List(MTSSLLabelingIndex), residue2List(MTSSLLabelingIndex), R1_20210523);
    end
    
    normalizedInitialSimulatedDistanceDistributionYPhaseThree = initialSimulatedDistanceDistributionYPhaseThree./sum(initialSimulatedDistanceDistributionYPhaseThree);
    
    [~, initialDEERMaximalIndexPhaseThree] = max(normalizedInitialSimulatedDistanceDistributionYPhaseThree);
    
    [~, targetDEERMaximalIndexTmp] = max(distanceDistributionList);
    
    targetDEERMaximalIndex = reshape(targetDEERMaximalIndexTmp(1, 2, :), 1, []);
    
    RMSEThirdPhase(1) = (mean((targetDEERMaximalIndex-initialDEERMaximalIndexPhaseThree).^2))^0.5;
    
    iterationIndexThirdPhase(1) = 1;

    monteCarloOldGeometryPhaseThree = monteCarloOldGeometryPhaseThree;
    monteCarloOldDistanceDistributionPhaseThree = initialDEERMaximalIndexPhaseThree;
    monteCarloOldStructurePhaseThree = monteCarloOldStructurePhaseThree;
    monteCarloOldContactedResidueIndexPhaseThree = formatedPDB2contactedResidueNumbers(monteCarloOldStructurePhaseThree, 2.4);
    clashedResidueNumberPhaseThree(1) = length(monteCarloOldContactedResidueIndexPhaseThree);
    for monteCarloStepsPhaseThree = 1:monteCarloIterationNumber
        fifthRMSEPassedNumber = load(fifthRMSEPassedNumberFileName);
        if fifthRMSEPassedNumber >= numberOfStructure
            return
        end
        monteCarloNewCandidateGeometryPhaseThree = monteCarloOldGeometryPhaseThree;
        monteCarloNewCandidateGeometryPhaseThree.dihedral.CN1CA1C1(flexibleRegionIndex) = monteCarloNewCandidateGeometryPhaseThree.dihedral.CN1CA1C1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoopPhaseThree;
        monteCarloNewCandidateGeometryPhaseThree.dihedral.NCACN1(flexibleRegionIndex) = monteCarloNewCandidateGeometryPhaseThree.dihedral.NCACN1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoopPhaseThree;
        currentCandidateBackboneCoordinatesPhaseThree = geometry2backboneCoordinates(monteCarloNewCandidateGeometryPhaseThree);
        [~, currentCandidateAlignedBackboneCoordinatesPhaseThree] = procrustes(initialBackboneCoordinates, currentCandidateBackboneCoordinatesPhaseThree, 'scaling', false);
        currentCandidateFormatedBackbonePhaseThree = initialFormatedBackbone;
        currentCandidateFormatedBackbonePhaseThree(:, 9:11) = string(currentCandidateAlignedBackboneCoordinatesPhaseThree);
%         currentCandidateMinimalBackboneNonBondedDistancePhaseThree = formatedBackbone2MinimalBackboneNonBondedDistance(currentCandidateFormatedBackbonePhaseThree);
%         if minimalBackboneNonBondedDistance*0.98 > currentCandidateMinimalBackboneNonBondedDistancePhaseThree 
%             currentStatusUpdator(currentStatusFileName, 6)
%             continue
%         end
        currentCandidateFormatedStructurePhaseThree = oldSideChainInstaller(currentCandidateFormatedBackbonePhaseThree, monteCarloOldStructurePhaseThree);
        newCandidateFormatedStructurePhaseThree_1 = sideChainRotator(currentCandidateFormatedStructurePhaseThree, 25, monteCarloOldContactedResidueIndexPhaseThree);
        newCandidateFormatedStructurePhaseThree = sideChainRotator(newCandidateFormatedStructurePhaseThree_1, 5, sideChainRotatorResidueIndex);
        newCandidateConteactedResidueIndexPhaseThree = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructurePhaseThree, 2.4);
        
        for labelingIndex = 1:length(residue1List)
            [newCandidateSimulatedDistanceDistributionYPhaseThree(:, labelingIndex), ~] = DEERefineMTSSLLabeling(newCandidateFormatedStructurePhaseThree, residue1List(labelingIndex), residue2List(labelingIndex), R1_20210523);
        end
        
        [~, newCandidateDEERMaximalIndexPhaseThree] = max(newCandidateSimulatedDistanceDistributionYPhaseThree);
        monteCarloNewDistanceDistributionPhaseThree = newCandidateDEERMaximalIndexPhaseThree;
        monteCarloClashesTemperaturePhaseThree = 0.92;
        [DEERIfMovingCriterion, RMSEThirdPhase(monteCarloStepsPhaseThree+1)] = monteCarloMetropolisCriterionGenerator(monteCarloOldDistanceDistributionPhaseThree, monteCarloNewDistanceDistributionPhaseThree, targetDEERMaximalIndex, monteCarloTemperature);
        [clashesIfMovingCriterion, clashedResidueNumberThirdPhase(monteCarloStepsPhaseThree+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndexPhaseThree), length(newCandidateConteactedResidueIndexPhaseThree), 0, monteCarloClashesTemperaturePhaseThree);
        
        if DEERIfMovingCriterion == 1 && clashesIfMovingCriterion == 1
            currentStatusUpdator(currentStatusFileName, 7)
            PDBEvolutionSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructurePhaseThree, PDBEvolution)
%             RMSEDatSaver(outputTrajectoryDatFileName, monteCarloStepsPhaseThree, RMSEThirdPhase(monteCarloStepsPhaseThree+1))
            datSaver(outputTrajectoryDatFileName, monteCarloStepsPhaseThree, RMSEThirdPhase(monteCarloStepsPhaseThree+1), clashedResidueNumberThirdPhase(monteCarloStepsPhaseThree+1))
            monteCarloOldContactedResidueIndexPhaseThree = newCandidateConteactedResidueIndexPhaseThree;
            monteCarloOldStructurePhaseThree = newCandidateFormatedStructurePhaseThree;
            monteCarloOldDistanceDistributionPhaseThree = monteCarloNewDistanceDistributionPhaseThree;
            monteCarloOldGeometryPhaseThree = monteCarloNewCandidateGeometryPhaseThree;
            if RMSEThirdPhase(monteCarloStepsPhaseThree+1) <= RMSE+0.2
                increaseSecondRMSEPassedNumber(secondRMSEPassedNumberFileName)

                currentStatusUpdator(currentStatusFileName, 8)
                monteCarloOldGeometryPhaseFour = monteCarloOldGeometryPhaseThree;
                monteCarloOldStructurePhaseFour = monteCarloOldStructurePhaseThree;
    %             pdbSaver(phaseOnePassedFileName, newCandidateFormatedStructure);
                break
            end
        else
            monteCarloOldContactedResidueIndexPhaseThree = monteCarloOldContactedResidueIndexPhaseThree;
            monteCarloOldStructurePhaseThree = monteCarloOldStructurePhaseThree;
            monteCarloOldDistanceDistributionPhaseThree = monteCarloOldDistanceDistributionPhaseThree;
            monteCarloOldGeometryPhaseThree = monteCarloOldGeometryPhaseThree;
        end
%         phaseOnePassedFileName = strcat(timeStamp, "_phase1_passed", num2str(simulationIndex), '_', num2str(structureIndex), '.pdb');

    end
    if exist('monteCarloOldStructurePhaseFour', 'var') ~= 1
        continue
    end
    clearvars -except simulationIndex structureIndex structureIndexEnd timeStamp...
        currentStatusFileName...
        firstRMSEPassedNumberFileName secondRMSEPassedNumberFileName thirdRMSEPassedNumberFileName fourthRMSEPassedNumberFileName fifthRMSEPassedNumberFileName...
        residue1List residue2List distanceDistributionList...
        initialBackboneCoordinates initialFormatedBackbone...
        sideChainRotatorResidueIndex minimalBackboneNonBondedDistance...
        phiPsiAngleVariationAtLoop flexibleRegionIndex...
        monteCarloOldGeometryPhaseFour monteCarloOldStructurePhaseFour...
        maximalClashes RMSE app numberOfStructure runFileName distanceDistributionFullPath...
                      R1LibraryMatFileName R1_20210523 PDBEvolution...
                      outputTrajectoryPDBFileName outputTrajectoryDatFileName monteCarloTemperature monteCarloIterationNumber
         
    
    
    
    
    phiPsiAngleVariationAtLoopPhaseFour = 0.05;
    clashesCriterionPhaseFour = 2.34;
    
    monteCarloOldContactedResidueIndexPhaseFour = formatedPDB2contactedResidueNumbers(monteCarloOldStructurePhaseFour, clashesCriterionPhaseFour);
    clashedResidueNumberPhaseFour(1) = length(monteCarloOldContactedResidueIndexPhaseFour);
    for phaseFourSteps = 1:monteCarloIterationNumber
        fifthRMSEPassedNumber = load(fifthRMSEPassedNumberFileName);
        if fifthRMSEPassedNumber >= numberOfStructure
            return
        end
        monteCarloNewCandidateGeometryPhaseFour = monteCarloOldGeometryPhaseFour;
        monteCarloNewCandidateGeometryPhaseFour.dihedral.CN1CA1C1(flexibleRegionIndex) = monteCarloNewCandidateGeometryPhaseFour.dihedral.CN1CA1C1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoopPhaseFour;
        monteCarloNewCandidateGeometryPhaseFour.dihedral.NCACN1(flexibleRegionIndex) = monteCarloNewCandidateGeometryPhaseFour.dihedral.NCACN1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoopPhaseFour;
        currentCandidateBackboneCoordinatesPhaseFour = geometry2backboneCoordinates(monteCarloNewCandidateGeometryPhaseFour);
        [~, currentCandidateAlignedBackboneCoordinatesPhaseFour] = procrustes(initialBackboneCoordinates, currentCandidateBackboneCoordinatesPhaseFour, 'scaling', false);
        currentCandidateFormatedBackbonePhaseFour = initialFormatedBackbone;
        currentCandidateFormatedBackbonePhaseFour(:, 9:11) = string(currentCandidateAlignedBackboneCoordinatesPhaseFour);
%         currentCandidateMinimalBackboneNonBondedDistance = formatedBackbone2MinimalBackboneNonBondedDistance(currentCandidateFormatedBackbonePhaseFour);
%         if minimalBackboneNonBondedDistance*0.98 > currentCandidateMinimalBackboneNonBondedDistance 
%             continue
%         end
        currentCandidateFormatedStructurePhaseFour = oldSideChainInstaller(currentCandidateFormatedBackbonePhaseFour, monteCarloOldStructurePhaseFour);
        newCandidateFormatedStructurePhaseFour = sideChainRotator(currentCandidateFormatedStructurePhaseFour, 25, monteCarloOldContactedResidueIndexPhaseFour);
        newCandidateConteactedResidueIndexPhaseFour = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructurePhaseFour, clashesCriterionPhaseFour);
        [clashesIfMovingCriterionPhaseFour, clashedResidueNumberPhaseFour(phaseFourSteps+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndexPhaseFour), length(newCandidateConteactedResidueIndexPhaseFour), 0, 0.6);
        if clashesIfMovingCriterionPhaseFour==1
            currentStatusUpdator(currentStatusFileName, 9)
            %%%%%%% print to status panel monte carlo moved clashes
            monteCarloOldContactedResidueIndexPhaseFour = newCandidateConteactedResidueIndexPhaseFour;
            monteCarloOldStructurePhaseFour = newCandidateFormatedStructurePhaseFour;
            monteCarloOldGeometryPhaseFour = monteCarloNewCandidateGeometryPhaseFour;
            PDBEvolutionSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructurePhaseFour, PDBEvolution)
            datSaver(outputTrajectoryDatFileName, phaseFourSteps, 0, clashedResidueNumberPhaseFour(phaseFourSteps+1))
            if clashedResidueNumberPhaseFour(phaseFourSteps+1) <= maximalClashes
                currentStatusUpdator(currentStatusFileName, 10)
                monteCarloOldGeometryPhaseFive = monteCarloNewCandidateGeometryPhaseFour;
                monteCarloOldStructurePhaseFive = newCandidateFormatedStructurePhaseFour;
                break
            end 
        else
            monteCarloOldContactedResidueIndexPhaseFour = monteCarloOldContactedResidueIndexPhaseFour;
            monteCarloOldStructurePhaseFour = monteCarloOldStructurePhaseFour;
            monteCarloOldGeometryPhaseFour = monteCarloOldGeometryPhaseFour;
        end
    end
    if exist('monteCarloOldStructurePhaseFive', 'var') ~= 1
        continue
    end
    
    
    
    clearvars -except simulationIndex structureIndex structureIndexEnd timeStamp...
        currentStatusFileName...
        firstRMSEPassedNumberFileName secondRMSEPassedNumberFileName thirdRMSEPassedNumberFileName fourthRMSEPassedNumberFileName fifthRMSEPassedNumberFileName...
        residue1List residue2List distanceDistributionList...
        initialBackboneCoordinates initialFormatedBackbone...
        sideChainRotatorResidueIndex minimalBackboneNonBondedDistance...
        phiPsiAngleVariationAtLoop flexibleRegionIndex...
        monteCarloOldGeometryPhaseFive monteCarloOldStructurePhaseFive...
        maximalClashes RMSE app numberOfStructure runFileName distanceDistributionFullPath...
                      R1LibraryMatFileName R1_20210523 PDBEvolution...
                      outputTrajectoryPDBFileName outputTrajectoryDatFileName monteCarloTemperature monteCarloIterationNumber
    phiPsiAngleVariationAtLoopPhaseFive = 0.05;
    
    for MTSSLLabelingIndex = 1:length(residue1List)
        [initialSimulatedDistanceDistributionYPhaseFive(:, MTSSLLabelingIndex), ~] = DEERefineMTSSLLabeling(monteCarloOldStructurePhaseFive, residue1List(MTSSLLabelingIndex), residue2List(MTSSLLabelingIndex), R1_20210523);
    end
    
    normalizedInitialSimulatedDistanceDistributionYPhaseFive = initialSimulatedDistanceDistributionYPhaseFive./sum(initialSimulatedDistanceDistributionYPhaseFive);
    
    [~, initialDEERMaximalIndexPhaseFive] = max(normalizedInitialSimulatedDistanceDistributionYPhaseFive);
    
    [~, targetDEERMaximalIndexTmp] = max(distanceDistributionList);
    
    targetDEERMaximalIndex = reshape(targetDEERMaximalIndexTmp(1, 2, :), 1, []);
    
    RMSEFifthPhase(1) = (mean((targetDEERMaximalIndex-initialDEERMaximalIndexPhaseFive).^2))^0.5;
    
    iterationIndexFifthPhase(1) = 1;

    monteCarloOldGeometryPhaseFive = monteCarloOldGeometryPhaseFive;
    monteCarloOldDistanceDistributionPhaseFive = initialDEERMaximalIndexPhaseFive;
    monteCarloOldStructurePhaseFive = monteCarloOldStructurePhaseFive;
    monteCarloOldContactedResidueIndexPhaseFive = formatedPDB2contactedResidueNumbers(monteCarloOldStructurePhaseFive, 2.4);
    clashedResidueNumberPhaseFive(1) = length(monteCarloOldContactedResidueIndexPhaseFive);
    
    for monteCarloStepsPhaseFive = 1:monteCarloIterationNumber
        fifthRMSEPassedNumber = load(fifthRMSEPassedNumberFileName);
        if fifthRMSEPassedNumber >= numberOfStructure
            return
        end
        monteCarloNewCandidateGeometryPhaseFive = monteCarloOldGeometryPhaseFive;
        monteCarloNewCandidateGeometryPhaseFive.dihedral.CN1CA1C1(flexibleRegionIndex) = monteCarloNewCandidateGeometryPhaseFive.dihedral.CN1CA1C1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoopPhaseFive;
        monteCarloNewCandidateGeometryPhaseFive.dihedral.NCACN1(flexibleRegionIndex) = monteCarloNewCandidateGeometryPhaseFive.dihedral.NCACN1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoopPhaseFive;
        currentCandidateBackboneCoordinatesPhaseFive = geometry2backboneCoordinates(monteCarloNewCandidateGeometryPhaseFive);
        [~, currentCandidateAlignedBackboneCoordinatesPhaseFive] = procrustes(initialBackboneCoordinates, currentCandidateBackboneCoordinatesPhaseFive, 'scaling', false);
        currentCandidateFormatedBackbonePhaseFive = initialFormatedBackbone;
        currentCandidateFormatedBackbonePhaseFive(:, 9:11) = string(currentCandidateAlignedBackboneCoordinatesPhaseFive);
%         currentCandidateMinimalBackboneNonBondedDistancePhaseFive = formatedBackbone2MinimalBackboneNonBondedDistance(currentCandidateFormatedBackbonePhaseFive);
%         if minimalBackboneNonBondedDistance*0.98 > currentCandidateMinimalBackboneNonBondedDistancePhaseFive 
%             currentStatusUpdator(currentStatusFileName, 11)
%             continue
%         end
        currentCandidateFormatedStructurePhaseFive = oldSideChainInstaller(currentCandidateFormatedBackbonePhaseFive, monteCarloOldStructurePhaseFive);
        newCandidateFormatedStructurePhaseFive_1 = sideChainRotator(currentCandidateFormatedStructurePhaseFive, 25, monteCarloOldContactedResidueIndexPhaseFive);
        newCandidateFormatedStructurePhaseFive = sideChainRotator(newCandidateFormatedStructurePhaseFive_1, 5, sideChainRotatorResidueIndex);
        newCandidateConteactedResidueIndexPhaseFive = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructurePhaseFive, 2.4);
        
        for labelingIndex = 1:length(residue1List)
            [newCandidateSimulatedDistanceDistributionYPhaseFive(:, labelingIndex), ~] = DEERefineMTSSLLabeling(newCandidateFormatedStructurePhaseFive, residue1List(labelingIndex), residue2List(labelingIndex), R1_20210523);
        end
        
        [~, newCandidateDEERMaximalIndexPhaseFive] = max(newCandidateSimulatedDistanceDistributionYPhaseFive);
        monteCarloNewDistanceDistributionPhaseFive = newCandidateDEERMaximalIndexPhaseFive;
        monteCarloClashesTemperaturePhaseFive = 0.92;
        [DEERIfMovingCriterion, RMSEFifthPhase(monteCarloStepsPhaseFive+1)] = monteCarloMetropolisCriterionGenerator(monteCarloOldDistanceDistributionPhaseFive, monteCarloNewDistanceDistributionPhaseFive, targetDEERMaximalIndex, monteCarloTemperature);
        [clashesIfMovingCriterion, clashedResidueNumberFifthPhase(monteCarloStepsPhaseFive+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndexPhaseFive), length(newCandidateConteactedResidueIndexPhaseFive), 0, monteCarloClashesTemperaturePhaseFive);
        
        if DEERIfMovingCriterion == 1 && clashesIfMovingCriterion == 1
            currentStatusUpdator(currentStatusFileName, 12)
            PDBEvolutionSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructurePhaseFive, PDBEvolution)
%             RMSEDatSaver(outputTrajectoryDatFileName, monteCarloStepsPhaseFive, RMSEFifthPhase(monteCarloStepsPhaseFive+1))
            datSaver(outputTrajectoryDatFileName, monteCarloStepsPhaseFive, RMSEFifthPhase(monteCarloStepsPhaseFive+1), clashedResidueNumberFifthPhase(monteCarloStepsPhaseFive+1))
            monteCarloOldContactedResidueIndexPhaseFive = newCandidateConteactedResidueIndexPhaseFive;
            monteCarloOldStructurePhaseFive = newCandidateFormatedStructurePhaseFive;
            monteCarloOldDistanceDistributionPhaseFive = monteCarloNewDistanceDistributionPhaseFive;
            monteCarloOldGeometryPhaseFive = monteCarloNewCandidateGeometryPhaseFive;
            if RMSEFifthPhase(monteCarloStepsPhaseFive+1) <= RMSE+0.2
                increaseThirdRMSEPassedNumber(thirdRMSEPassedNumberFileName)

                currentStatusUpdator(currentStatusFileName, 13)
                monteCarloOldGeometryPhaseSix = monteCarloOldGeometryPhaseFive;
                monteCarloOldStructurePhaseSix = monteCarloOldStructurePhaseFive;
    %             pdbSaver(phaseOnePassedFileName, newCandidateFormatedStructure);
                break
            end 
        else
            monteCarloOldContactedResidueIndexPhaseFive = monteCarloOldContactedResidueIndexPhaseFive;
            monteCarloOldStructurePhaseFive = monteCarloOldStructurePhaseFive;
            monteCarloOldDistanceDistributionPhaseFive = monteCarloOldDistanceDistributionPhaseFive;
            monteCarloOldGeometryPhaseFive = monteCarloOldGeometryPhaseFive;
        end
%         phaseOnePassedFileName = strcat(timeStamp, "_phase1_passed", num2str(simulationIndex), '_', num2str(structureIndex), '.pdb');

    end
    if exist('monteCarloOldStructurePhaseSix', 'var') ~= 1
        continue
    end
    
    
    
    
    clearvars -except simulationIndex structureIndex structureIndexEnd timeStamp...
        currentStatusFileName...
        firstRMSEPassedNumberFileName secondRMSEPassedNumberFileName thirdRMSEPassedNumberFileName fourthRMSEPassedNumberFileName fifthRMSEPassedNumberFileName...
        residue1List residue2List distanceDistributionList...
        initialBackboneCoordinates initialFormatedBackbone...
        sideChainRotatorResidueIndex minimalBackboneNonBondedDistance...
        phiPsiAngleVariationAtLoop flexibleRegionIndex...
        monteCarloOldGeometryPhaseSix monteCarloOldStructurePhaseSix...
        maximalClashes RMSE app numberOfStructure runFileName distanceDistributionFullPath...
                      R1LibraryMatFileName R1_20210523 PDBEvolution...
                      outputTrajectoryPDBFileName outputTrajectoryDatFileName monteCarloTemperature monteCarloIterationNumber
    
    
    
    
    phiPsiAngleVariationAtLoopPhaseSix = 0.05;
    clashesCriterionPhaseSix = 2.38;
    
    monteCarloOldContactedResidueIndexPhaseSix = formatedPDB2contactedResidueNumbers(monteCarloOldStructurePhaseSix, clashesCriterionPhaseSix);
    clashedResidueNumberPhaseSix(1) = length(monteCarloOldContactedResidueIndexPhaseSix);
    for phaseSixSteps = 1:monteCarloIterationNumber
        fifthRMSEPassedNumber = load(fifthRMSEPassedNumberFileName);
        if fifthRMSEPassedNumber >= numberOfStructure
            return
        end
        monteCarloNewCandidateGeometryPhaseSix = monteCarloOldGeometryPhaseSix;
        monteCarloNewCandidateGeometryPhaseSix.dihedral.CN1CA1C1(flexibleRegionIndex) = monteCarloNewCandidateGeometryPhaseSix.dihedral.CN1CA1C1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoopPhaseSix;
        monteCarloNewCandidateGeometryPhaseSix.dihedral.NCACN1(flexibleRegionIndex) = monteCarloNewCandidateGeometryPhaseSix.dihedral.NCACN1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoopPhaseSix;
        currentCandidateBackboneCoordinatesPhaseSix = geometry2backboneCoordinates(monteCarloNewCandidateGeometryPhaseSix);
        [~, currentCandidateAlignedBackboneCoordinatesPhaseSix] = procrustes(initialBackboneCoordinates, currentCandidateBackboneCoordinatesPhaseSix, 'scaling', false);
        currentCandidateFormatedBackbonePhaseSix = initialFormatedBackbone;
        currentCandidateFormatedBackbonePhaseSix(:, 9:11) = string(currentCandidateAlignedBackboneCoordinatesPhaseSix);
%         currentCandidateMinimalBackboneNonBondedDistance = formatedBackbone2MinimalBackboneNonBondedDistance(currentCandidateFormatedBackbonePhaseSix);
%         if minimalBackboneNonBondedDistance*0.98 > currentCandidateMinimalBackboneNonBondedDistance 
%             continue
%         end
        currentCandidateFormatedStructurePhaseSix = oldSideChainInstaller(currentCandidateFormatedBackbonePhaseSix, monteCarloOldStructurePhaseSix);
        newCandidateFormatedStructurePhaseSix = sideChainRotator(currentCandidateFormatedStructurePhaseSix, 25, monteCarloOldContactedResidueIndexPhaseSix);
        newCandidateConteactedResidueIndexPhaseSix = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructurePhaseSix, clashesCriterionPhaseSix);
        [clashesIfMovingCriterionPhaseSix, clashedResidueNumberPhaseSix(phaseSixSteps+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndexPhaseSix), length(newCandidateConteactedResidueIndexPhaseSix), 0, 0.6);
        if clashesIfMovingCriterionPhaseSix==1
            currentStatusUpdator(currentStatusFileName, 14)
            monteCarloOldContactedResidueIndexPhaseSix = newCandidateConteactedResidueIndexPhaseSix;
            monteCarloOldStructurePhaseSix = newCandidateFormatedStructurePhaseSix;
            monteCarloOldGeometryPhaseSix = monteCarloNewCandidateGeometryPhaseSix;
            PDBEvolutionSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructurePhaseSix, PDBEvolution)
            datSaver(outputTrajectoryDatFileName, phaseSixSteps, 0, clashedResidueNumberPhaseSix(phaseSixSteps+1))
            if clashedResidueNumberPhaseSix(phaseSixSteps+1) <= maximalClashes
                currentStatusUpdator(currentStatusFileName, 15)
                monteCarloOldGeometryPhaseSeven = monteCarloNewCandidateGeometryPhaseSix;
                monteCarloOldStructurePhaseSeven = newCandidateFormatedStructurePhaseSix;
                break
            end 
        else
            monteCarloOldContactedResidueIndexPhaseSix = monteCarloOldContactedResidueIndexPhaseSix;
            monteCarloOldStructurePhaseSix = monteCarloOldStructurePhaseSix;
            monteCarloOldGeometryPhaseSix = monteCarloOldGeometryPhaseSix;
        end
    end
    if exist('monteCarloOldStructurePhaseSeven', 'var') ~= 1
        continue
    end
    
    
    
    
    
    clearvars -except simulationIndex structureIndex structureIndexEnd timeStamp...
        currentStatusFileName...
        firstRMSEPassedNumberFileName secondRMSEPassedNumberFileName thirdRMSEPassedNumberFileName fourthRMSEPassedNumberFileName fifthRMSEPassedNumberFileName...
        residue1List residue2List distanceDistributionList...
        initialBackboneCoordinates initialFormatedBackbone...
        sideChainRotatorResidueIndex minimalBackboneNonBondedDistance...
        phiPsiAngleVariationAtLoop flexibleRegionIndex...
        monteCarloOldGeometryPhaseSeven monteCarloOldStructurePhaseSeven...
        maximalClashes RMSE app numberOfStructure runFileName distanceDistributionFullPath...
                      R1LibraryMatFileName R1_20210523 PDBEvolution...
                      outputTrajectoryPDBFileName outputTrajectoryDatFileName monteCarloTemperature monteCarloIterationNumber
    phiPsiAngleVariationAtLoopPhaseSeven = 0.05;
    
    for MTSSLLabelingIndex = 1:length(residue1List)
        [initialSimulatedDistanceDistributionYPhaseSeven(:, MTSSLLabelingIndex), ~] = DEERefineMTSSLLabeling(monteCarloOldStructurePhaseSeven, residue1List(MTSSLLabelingIndex), residue2List(MTSSLLabelingIndex), R1_20210523);
    end
    
    normalizedInitialSimulatedDistanceDistributionYPhaseSeven = initialSimulatedDistanceDistributionYPhaseSeven./sum(initialSimulatedDistanceDistributionYPhaseSeven);
    
    [~, initialDEERMaximalIndexPhaseSeven] = max(normalizedInitialSimulatedDistanceDistributionYPhaseSeven);
    
    [~, targetDEERMaximalIndexTmp] = max(distanceDistributionList);
    
    targetDEERMaximalIndex = reshape(targetDEERMaximalIndexTmp(1, 2, :), 1, []);
    
    RMSESeventhPhase(1) = (mean((targetDEERMaximalIndex-initialDEERMaximalIndexPhaseSeven).^2))^0.5;
    
    iterationIndexSeventhPhase(1) = 1;

    monteCarloOldGeometryPhaseSeven = monteCarloOldGeometryPhaseSeven;
    monteCarloOldDistanceDistributionPhaseSeven = initialDEERMaximalIndexPhaseSeven;
    monteCarloOldStructurePhaseSeven = monteCarloOldStructurePhaseSeven;
    monteCarloOldContactedResidueIndexPhaseSeven = formatedPDB2contactedResidueNumbers(monteCarloOldStructurePhaseSeven, 2.4);
    clashedResidueNumberPhaseSeven(1) = length(monteCarloOldContactedResidueIndexPhaseSeven);
    
    for monteCarloStepsPhaseSeven = 1:monteCarloIterationNumber
        fifthRMSEPassedNumber = load(fifthRMSEPassedNumberFileName);
        if fifthRMSEPassedNumber >= numberOfStructure
            return
        end
        monteCarloNewCandidateGeometryPhaseSeven = monteCarloOldGeometryPhaseSeven;
        monteCarloNewCandidateGeometryPhaseSeven.dihedral.CN1CA1C1(flexibleRegionIndex) = monteCarloNewCandidateGeometryPhaseSeven.dihedral.CN1CA1C1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoopPhaseSeven;
        monteCarloNewCandidateGeometryPhaseSeven.dihedral.NCACN1(flexibleRegionIndex) = monteCarloNewCandidateGeometryPhaseSeven.dihedral.NCACN1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoopPhaseSeven;
        currentCandidateBackboneCoordinatesPhaseSeven = geometry2backboneCoordinates(monteCarloNewCandidateGeometryPhaseSeven);
        [~, currentCandidateAlignedBackboneCoordinatesPhaseSeven] = procrustes(initialBackboneCoordinates, currentCandidateBackboneCoordinatesPhaseSeven, 'scaling', false);
        currentCandidateFormatedBackbonePhaseSeven = initialFormatedBackbone;
        currentCandidateFormatedBackbonePhaseSeven(:, 9:11) = string(currentCandidateAlignedBackboneCoordinatesPhaseSeven);
%         currentCandidateMinimalBackboneNonBondedDistancePhaseSeven = formatedBackbone2MinimalBackboneNonBondedDistance(currentCandidateFormatedBackbonePhaseSeven);
%         if minimalBackboneNonBondedDistance*0.98 > currentCandidateMinimalBackboneNonBondedDistancePhaseSeven 
%             currentStatusUpdator(currentStatusFileName, 16)
%             continue
%         end
        currentCandidateFormatedStructurePhaseSeven = oldSideChainInstaller(currentCandidateFormatedBackbonePhaseSeven, monteCarloOldStructurePhaseSeven);
        newCandidateFormatedStructurePhaseSeven_1 = sideChainRotator(currentCandidateFormatedStructurePhaseSeven, 25, monteCarloOldContactedResidueIndexPhaseSeven);
        newCandidateFormatedStructurePhaseSeven = sideChainRotator(newCandidateFormatedStructurePhaseSeven_1, 5, sideChainRotatorResidueIndex);
        newCandidateConteactedResidueIndexPhaseSeven = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructurePhaseSeven, 2.4);
        
        for labelingIndex = 1:length(residue1List)
            [newCandidateSimulatedDistanceDistributionYPhaseSeven(:, labelingIndex), ~] = DEERefineMTSSLLabeling(newCandidateFormatedStructurePhaseSeven, residue1List(labelingIndex), residue2List(labelingIndex), R1_20210523);
        end
        
        [~, newCandidateDEERMaximalIndexPhaseSeven] = max(newCandidateSimulatedDistanceDistributionYPhaseSeven);
        monteCarloNewDistanceDistributionPhaseSeven = newCandidateDEERMaximalIndexPhaseSeven;
        monteCarloClashesTemperaturePhaseSeven = 0.92;
        [DEERIfMovingCriterion, RMSESeventhPhase(monteCarloStepsPhaseSeven+1)] = monteCarloMetropolisCriterionGenerator(monteCarloOldDistanceDistributionPhaseSeven, monteCarloNewDistanceDistributionPhaseSeven, targetDEERMaximalIndex, monteCarloTemperature);
        [clashesIfMovingCriterion, clashedResidueNumberSeventhPhase(monteCarloStepsPhaseSeven+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndexPhaseSeven), length(newCandidateConteactedResidueIndexPhaseSeven), 0, monteCarloClashesTemperaturePhaseSeven);
        
        if DEERIfMovingCriterion == 1 && clashesIfMovingCriterion == 1
            currentStatusUpdator(currentStatusFileName, 17)
            PDBEvolutionSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructurePhaseSeven, PDBEvolution)
%             RMSEDatSaver(outputTrajectoryDatFileName, monteCarloStepsPhaseSeven, RMSESeventhPhase(monteCarloStepsPhaseSeven+1))
            datSaver(outputTrajectoryDatFileName, monteCarloStepsPhaseSeven, RMSESeventhPhase(monteCarloStepsPhaseSeven+1), clashedResidueNumberSeventhPhase(monteCarloStepsPhaseSeven+1))
            monteCarloOldContactedResidueIndexPhaseSeven = newCandidateConteactedResidueIndexPhaseSeven;
            monteCarloOldStructurePhaseSeven = newCandidateFormatedStructurePhaseSeven;
            monteCarloOldDistanceDistributionPhaseSeven = monteCarloNewDistanceDistributionPhaseSeven;
            monteCarloOldGeometryPhaseSeven = monteCarloNewCandidateGeometryPhaseSeven;
            if RMSESeventhPhase(monteCarloStepsPhaseSeven+1) <= RMSE+0.2
                increaseFourthRMSEPassedNumber(fourthRMSEPassedNumberFileName)

                currentStatusUpdator(currentStatusFileName, 18)
                monteCarloOldGeometryPhaseEight = monteCarloOldGeometryPhaseSeven;
                monteCarloOldStructurePhaseEight = monteCarloOldStructurePhaseSeven;
    %             pdbSaver(phaseOnePassedFileName, newCandidateFormatedStructure);
                break
            end 
        else
            monteCarloOldContactedResidueIndexPhaseSeven = monteCarloOldContactedResidueIndexPhaseSeven;
            monteCarloOldStructurePhaseSeven = monteCarloOldStructurePhaseSeven;
            monteCarloOldDistanceDistributionPhaseSeven = monteCarloOldDistanceDistributionPhaseSeven;
            monteCarloOldGeometryPhaseSeven = monteCarloOldGeometryPhaseSeven;
        end
%         phaseOnePassedFileName = strcat(timeStamp, "_phase1_passed", num2str(simulationIndex), '_', num2str(structureIndex), '.pdb');

    end
    if exist('monteCarloOldStructurePhaseEight', 'var') ~= 1
        continue
    end
    
    clearvars -except simulationIndex structureIndex structureIndexEnd timeStamp...
        currentStatusFileName...
        firstRMSEPassedNumberFileName secondRMSEPassedNumberFileName thirdRMSEPassedNumberFileName fourthRMSEPassedNumberFileName fifthRMSEPassedNumberFileName...
        residue1List residue2List distanceDistributionList...
        initialBackboneCoordinates initialFormatedBackbone...
        sideChainRotatorResidueIndex minimalBackboneNonBondedDistance...
        phiPsiAngleVariationAtLoop flexibleRegionIndex...
        monteCarloOldGeometryPhaseEight monteCarloOldStructurePhaseEight...
        maximalClashes RMSE app numberOfStructure runFileName distanceDistributionFullPath...
                      R1LibraryMatFileName R1_20210523 PDBEvolution...
                      outputTrajectoryPDBFileName outputTrajectoryDatFileName monteCarloTemperature monteCarloIterationNumber
    
    
    
    
    phiPsiAngleVariationAtLoopPhaseEight = 0.05;
    clashesCriterionPhaseEight = 2.42;
    
    monteCarloOldContactedResidueIndexPhaseEight = formatedPDB2contactedResidueNumbers(monteCarloOldStructurePhaseEight, clashesCriterionPhaseEight);
    clashedResidueNumberPhaseEight(1) = length(monteCarloOldContactedResidueIndexPhaseEight);
    for phaseEightSteps = 1:monteCarloIterationNumber
        fifthRMSEPassedNumber = load(fifthRMSEPassedNumberFileName);
        if fifthRMSEPassedNumber >= numberOfStructure
            return
        end
        monteCarloNewCandidateGeometryPhaseEight = monteCarloOldGeometryPhaseEight;
        monteCarloNewCandidateGeometryPhaseEight.dihedral.CN1CA1C1(flexibleRegionIndex) = monteCarloNewCandidateGeometryPhaseEight.dihedral.CN1CA1C1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoopPhaseEight;
        monteCarloNewCandidateGeometryPhaseEight.dihedral.NCACN1(flexibleRegionIndex) = monteCarloNewCandidateGeometryPhaseEight.dihedral.NCACN1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoopPhaseEight;
        currentCandidateBackboneCoordinatesPhaseEight = geometry2backboneCoordinates(monteCarloNewCandidateGeometryPhaseEight);
        [~, currentCandidateAlignedBackboneCoordinatesPhaseEight] = procrustes(initialBackboneCoordinates, currentCandidateBackboneCoordinatesPhaseEight, 'scaling', false);
        currentCandidateFormatedBackbonePhaseEight = initialFormatedBackbone;
        currentCandidateFormatedBackbonePhaseEight(:, 9:11) = string(currentCandidateAlignedBackboneCoordinatesPhaseEight);
%         currentCandidateMinimalBackboneNonBondedDistance = formatedBackbone2MinimalBackboneNonBondedDistance(currentCandidateFormatedBackbonePhaseEight);
%         if minimalBackboneNonBondedDistance*0.98 > currentCandidateMinimalBackboneNonBondedDistance 
%             continue
%         end
        currentCandidateFormatedStructurePhaseEight = oldSideChainInstaller(currentCandidateFormatedBackbonePhaseEight, monteCarloOldStructurePhaseEight);
        newCandidateFormatedStructurePhaseEight = sideChainRotator(currentCandidateFormatedStructurePhaseEight, 25, monteCarloOldContactedResidueIndexPhaseEight);
        newCandidateConteactedResidueIndexPhaseEight = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructurePhaseEight, clashesCriterionPhaseEight);
        [clashesIfMovingCriterionPhaseEight, clashedResidueNumberPhaseEight(phaseEightSteps+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndexPhaseEight), length(newCandidateConteactedResidueIndexPhaseEight), 0, 0.6);
        if clashesIfMovingCriterionPhaseEight==1
            currentStatusUpdator(currentStatusFileName, 19)
            monteCarloOldContactedResidueIndexPhaseEight = newCandidateConteactedResidueIndexPhaseEight;
            monteCarloOldStructurePhaseEight = newCandidateFormatedStructurePhaseEight;
            monteCarloOldGeometryPhaseEight = monteCarloNewCandidateGeometryPhaseEight;
            PDBEvolutionSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructurePhaseEight, PDBEvolution)
            datSaver(outputTrajectoryDatFileName, phaseEightSteps, 0, clashedResidueNumberPhaseEight(phaseEightSteps+1))
            if clashedResidueNumberPhaseEight(phaseEightSteps+1) <= maximalClashes
                currentStatusUpdator(currentStatusFileName, 20)
                monteCarloOldGeometryPhaseNine = monteCarloNewCandidateGeometryPhaseEight;
                monteCarloOldStructurePhaseNine = newCandidateFormatedStructurePhaseEight;
                break
            end 
        else
            monteCarloOldContactedResidueIndexPhaseEight = monteCarloOldContactedResidueIndexPhaseEight;
            monteCarloOldStructurePhaseEight = monteCarloOldStructurePhaseEight;
            monteCarloOldGeometryPhaseEight = monteCarloOldGeometryPhaseEight;
        end
    end
    if exist('monteCarloOldStructurePhaseNine', 'var') ~= 1
        continue
    end
    %%%%%
    clearvars -except simulationIndex structureIndex structureIndexEnd timeStamp...
        currentStatusFileName...
        firstRMSEPassedNumberFileName secondRMSEPassedNumberFileName thirdRMSEPassedNumberFileName fourthRMSEPassedNumberFileName fifthRMSEPassedNumberFileName...
        residue1List residue2List distanceDistributionList...
        initialBackboneCoordinates initialFormatedBackbone...
        sideChainRotatorResidueIndex minimalBackboneNonBondedDistance...
        phiPsiAngleVariationAtLoop flexibleRegionIndex...
        monteCarloOldGeometryPhaseNine monteCarloOldStructurePhaseNine...
        maximalClashes RMSE app numberOfStructure runFileName distanceDistributionFullPath...
                      R1LibraryMatFileName R1_20210523 PDBEvolution...
                      outputTrajectoryPDBFileName outputTrajectoryDatFileName monteCarloTemperature monteCarloIterationNumber
    phiPsiAngleVariationAtLoopPhaseNine = 0.05;
    
    for MTSSLLabelingIndex = 1:length(residue1List)
        [initialSimulatedDistanceDistributionYPhaseNine(:, MTSSLLabelingIndex), ~] = DEERefineMTSSLLabeling(monteCarloOldStructurePhaseNine, residue1List(MTSSLLabelingIndex), residue2List(MTSSLLabelingIndex), R1_20210523);
    end
    
    normalizedInitialSimulatedDistanceDistributionYPhaseNine = initialSimulatedDistanceDistributionYPhaseNine./sum(initialSimulatedDistanceDistributionYPhaseNine);
    
    [~, initialDEERMaximalIndexPhaseNine] = max(normalizedInitialSimulatedDistanceDistributionYPhaseNine);
    
    [~, targetDEERMaximalIndexTmp] = max(distanceDistributionList);
    
    targetDEERMaximalIndex = reshape(targetDEERMaximalIndexTmp(1, 2, :), 1, []);
    
    RMSENinethPhase(1) = (mean((targetDEERMaximalIndex-initialDEERMaximalIndexPhaseNine).^2))^0.5;
    
    iterationIndexNinethPhase(1) = 1;

    monteCarloOldGeometryPhaseNine = monteCarloOldGeometryPhaseNine;
    monteCarloOldDistanceDistributionPhaseNine = initialDEERMaximalIndexPhaseNine;
    monteCarloOldStructurePhaseNine = monteCarloOldStructurePhaseNine;
    monteCarloOldContactedResidueIndexPhaseNine = formatedPDB2contactedResidueNumbers(monteCarloOldStructurePhaseNine, 2.4);
    clashedResidueNumberPhaseNine(1) = length(monteCarloOldContactedResidueIndexPhaseNine);
    
    for monteCarloStepsPhaseNine = 1:monteCarloIterationNumber
        fifthRMSEPassedNumber = load(fifthRMSEPassedNumberFileName);
        if fifthRMSEPassedNumber >= numberOfStructure
            return
        end
        monteCarloNewCandidateGeometryPhaseNine = monteCarloOldGeometryPhaseNine;
        monteCarloNewCandidateGeometryPhaseNine.dihedral.CN1CA1C1(flexibleRegionIndex) = monteCarloNewCandidateGeometryPhaseNine.dihedral.CN1CA1C1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoopPhaseNine;
        monteCarloNewCandidateGeometryPhaseNine.dihedral.NCACN1(flexibleRegionIndex) = monteCarloNewCandidateGeometryPhaseNine.dihedral.NCACN1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoopPhaseNine;
        currentCandidateBackboneCoordinatesPhaseNine = geometry2backboneCoordinates(monteCarloNewCandidateGeometryPhaseNine);
        [~, currentCandidateAlignedBackboneCoordinatesPhaseNine] = procrustes(initialBackboneCoordinates, currentCandidateBackboneCoordinatesPhaseNine, 'scaling', false);
        currentCandidateFormatedBackbonePhaseNine = initialFormatedBackbone;
        currentCandidateFormatedBackbonePhaseNine(:, 9:11) = string(currentCandidateAlignedBackboneCoordinatesPhaseNine);
%         currentCandidateMinimalBackboneNonBondedDistancePhaseNine = formatedBackbone2MinimalBackboneNonBondedDistance(currentCandidateFormatedBackbonePhaseNine);
%         if minimalBackboneNonBondedDistance*0.98 > currentCandidateMinimalBackboneNonBondedDistancePhaseNine 
%             currentStatusUpdator(currentStatusFileName, 21)
%             continue
%         end
        currentCandidateFormatedStructurePhaseNine = oldSideChainInstaller(currentCandidateFormatedBackbonePhaseNine, monteCarloOldStructurePhaseNine);
        newCandidateFormatedStructurePhaseNine_1 = sideChainRotator(currentCandidateFormatedStructurePhaseNine, 25, monteCarloOldContactedResidueIndexPhaseNine);
        newCandidateFormatedStructurePhaseNine = sideChainRotator(newCandidateFormatedStructurePhaseNine_1, 5, sideChainRotatorResidueIndex);
        newCandidateConteactedResidueIndexPhaseNine = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructurePhaseNine, 2.4);
        
        for labelingIndex = 1:length(residue1List)
            [newCandidateSimulatedDistanceDistributionYPhaseNine(:, labelingIndex), ~] = DEERefineMTSSLLabeling(newCandidateFormatedStructurePhaseNine, residue1List(labelingIndex), residue2List(labelingIndex), R1_20210523);
        end
        
        [~, newCandidateDEERMaximalIndexPhaseNine] = max(newCandidateSimulatedDistanceDistributionYPhaseNine);
        monteCarloNewDistanceDistributionPhaseNine = newCandidateDEERMaximalIndexPhaseNine;
        monteCarloClashesTemperaturePhaseNine = 0.85;
        [DEERIfMovingCriterion, RMSENinethPhase(monteCarloStepsPhaseNine+1)] = monteCarloMetropolisCriterionGenerator(monteCarloOldDistanceDistributionPhaseNine, monteCarloNewDistanceDistributionPhaseNine, targetDEERMaximalIndex, monteCarloTemperature);
        [clashesIfMovingCriterion, clashedResidueNumberNinethPhase(monteCarloStepsPhaseNine+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndexPhaseNine), length(newCandidateConteactedResidueIndexPhaseNine), 0, monteCarloClashesTemperaturePhaseNine);
        
        if DEERIfMovingCriterion == 1 && clashesIfMovingCriterion == 1
            currentStatusUpdator(currentStatusFileName, 22)
            PDBEvolutionSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructurePhaseNine, PDBEvolution)
%             RMSEDatSaver(outputTrajectoryDatFileName, monteCarloStepsPhaseNine, RMSENinethPhase(monteCarloStepsPhaseNine+1))
            datSaver(outputTrajectoryDatFileName, monteCarloStepsPhaseNine, RMSENinethPhase(monteCarloStepsPhaseNine+1), clashedResidueNumberNinethPhase(monteCarloStepsPhaseNine+1))
            monteCarloOldContactedResidueIndexPhaseNine = newCandidateConteactedResidueIndexPhaseNine;
            monteCarloOldStructurePhaseNine = newCandidateFormatedStructurePhaseNine;
            monteCarloOldDistanceDistributionPhaseNine = monteCarloNewDistanceDistributionPhaseNine;
            monteCarloOldGeometryPhaseNine = monteCarloNewCandidateGeometryPhaseNine;
            maximalAbsoluteRMSEDifference = max(abs(monteCarloOldDistanceDistributionPhaseNine-targetDEERMaximalIndex));
            if maximalAbsoluteRMSEDifference <= round(RMSE*2) && clashedResidueNumberNinethPhase(monteCarloStepsPhaseNine+1) <= maximalClashes*1.8 && RMSENinethPhase(monteCarloStepsPhaseNine+1) <= RMSE+0.4
                increaseFifthRMSEPassedNumber(fifthRMSEPassedNumberFileName)

                fifthRMSEPassedNumber = load(fifthRMSEPassedNumberFileName);
                candidateGenerator(fifthRMSEPassedNumber, newCandidateFormatedStructurePhaseNine, timeStamp)

                currentStatusUpdator(currentStatusFileName, 23)

                if fifthRMSEPassedNumber == numberOfStructure
                    FINALPDBENSEMBLEFileGenerator(timeStamp)
                    pause(5)
                    FINALPDBENSEMBLEFile = dir('*FINALPDBENSEMBLE_*');
                    FINALPDBENSEMBLEFileName = FINALPDBENSEMBLEFile.name;


                    [distanceDistributionFinal, FINALPDBENSEMBLE] = distanceDistributionFinalGenerator(FINALPDBENSEMBLEFileName, residue1List, residue2List, R1_20210523);

                    [~, targetDistanceDistributionFinal] = targetDistanceDistributionFinalForJSDGenerator(distanceDistributionFullPath);



                    
                    for FINALPDBNumber = 1:length(FINALPDBENSEMBLE(1, 1, :))
                        [~, FINALPDBMaximalIndex(FINALPDBNumber, :)] = max(distanceDistributionFinal(:, :, FINALPDBNumber));
                    end
                    differenceBetweenFINALPDBAndTarget = abs(FINALPDBMaximalIndex-targetDEERMaximalIndex);
                    FINALPDBMaximalDifference = max(differenceBetweenFINALPDBAndTarget');
                    FINALPDBRMSE = sqrt(mean((differenceBetweenFINALPDBAndTarget.^2)'));
                    
                    minimalFINALPDBMaximalDifference = min(FINALPDBMaximalDifference);
                    minimalFINALPDBRMSE = min(FINALPDBRMSE(FINALPDBMaximalDifference == min(FINALPDBMaximalDifference)));
                    

                    minimalRMSEIndex = minimalRMSEPDBGenerator(FINALPDBENSEMBLE, FINALPDBMaximalDifference, minimalFINALPDBMaximalDifference, FINALPDBRMSE, minimalFINALPDBRMSE, timeStamp);

                    DEERefineFinalFileGenerator(residue1List, residue2List, FINALPDBENSEMBLE, minimalRMSEIndex, targetDistanceDistributionFinal, distanceDistributionFullPath, timeStamp, R1_20210523)

                    DEERefineFinalFileCopier(timeStamp)
                    closereq
                    return
                elseif fifthRMSEPassedNumber > numberOfStructure
                    return
                end

                break
            end
        else
            monteCarloOldContactedResidueIndexPhaseNine = monteCarloOldContactedResidueIndexPhaseNine;
            monteCarloOldStructurePhaseNine = monteCarloOldStructurePhaseNine;
            monteCarloOldDistanceDistributionPhaseNine = monteCarloOldDistanceDistributionPhaseNine;
            monteCarloOldGeometryPhaseNine = monteCarloOldGeometryPhaseNine;
        end
%         phaseOnePassedFileName = strcat(timeStamp, "_phase1_passed", num2str(simulationIndex), '_', num2str(structureIndex), '.pdb');

    end

end
