R1LibraryMatFileName = 'R1_library_20210523.mat';
load(R1LibraryMatFileName);
DEERefineSturctureFileName = "DEERefineSturcture.mat";
load(DEERefineSturctureFileName);
structureIndexEnd = structureIndexEachScript;
for structureIndex = 1:structureIndexEnd
    clearvars -except structureIndex simulationIndex structureIndexEnd app ...
                      runFileName RMSE maximalClashes monteCarloIteration initialStructureFullPath phiPsiAngle flexiblePhiPsiIndex targetNumberStructure...
                      R1LibraryMatFileName R1_20210523
                  
    DEERefineSturctureFileName = "DEERefineSturcture.mat";
    load(DEERefineSturctureFileName);
    rng('shuffle')
    timeStamp = char(erase(runFileName(1),"runFile_"));
    timeStamp(end-1:end) = [];
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
    
    outputTrajectoryPDBFileName = strcat(timeStamp, '_Trajectory_', num2str(1), '_', num2str(structureIndex), '.pdb');
    
    normalizedInitialSimulatedDistanceDistributionY = initialSimulatedDistanceDistributionY./sum(initialSimulatedDistanceDistributionY);
    
    [~, targetDEERMaximalIndexTmp] = max(distanceDistributionList);
    
    targetDEERMaximalIndex = reshape(targetDEERMaximalIndexTmp(1, 2, :), 1, []);
    
    [~, initialDEERMaximalIndex] = max(normalizedInitialSimulatedDistanceDistributionY);
    
    RMSEFirstPhase(1) = (mean((targetDEERMaximalIndex-initialDEERMaximalIndex).^2))^0.5;
    
    iterationIndexFirstPhase(1) = 1;
    
    outputTrajectoryDatFileName = strcat(timeStamp, '_Trajectory_', num2str(1), '_', num2str(structureIndex), '.dat');
    
    phiPsiAngleVariationAtLoop = phiPsiAngle;
    
    allResidueIndex = unique(double(initialMutatedFormatedPDB(:, 7)));

    flexibleRegionIndex = flexiblePhiPsiIndex;
    
    monteCarloOldGeometry = initialBackboneGeometry;
    monteCarloOldDistanceDistribution = initialDEERMaximalIndex;
    monteCarloOldStructure = initialMutatedFormatedPDB;
    monteCarloOldContactedResidueIndex = formatedPDB2contactedResidueNumbers(monteCarloOldStructure, 2.351);
    
    clashedResidueNumber(1) = length(monteCarloOldContactedResidueIndex);
    
    for monteCarloSteps = 1:monteCarloIterationNumber
        monteCarloNewCandidateGeometry = monteCarloOldGeometry;
        monteCarloNewCandidateGeometry.dihedral.CN1CA1C1(flexibleRegionIndex) = monteCarloNewCandidateGeometry.dihedral.CN1CA1C1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoop;
        monteCarloNewCandidateGeometry.dihedral.NCACN1(flexibleRegionIndex) = monteCarloNewCandidateGeometry.dihedral.NCACN1(flexibleRegionIndex)+(rand(length(flexibleRegionIndex), 1)-0.5)*phiPsiAngleVariationAtLoop;
        currentCandidateBackboneCoordinates = geometry2backboneCoordinates(monteCarloNewCandidateGeometry);
        [~, currentCandidateAlignedBackboneCoordinates] = procrustes(initialBackboneCoordinates, currentCandidateBackboneCoordinates, 'scaling', false);
        currentCandidateFormatedBackbone = initialFormatedBackbone;
        currentCandidateFormatedBackbone(:, 9:11) = string(currentCandidateAlignedBackboneCoordinates);
        currentCandidateMinimalBackboneNonBondedDistance = formatedBackbone2MinimalBackboneNonBondedDistance(currentCandidateFormatedBackbone);
        if minimalBackboneNonBondedDistance*0.98 > currentCandidateMinimalBackboneNonBondedDistance 
            currentStatusUpdator(currentStatusFileName, 1)
            continue
        end
        currentCandidateFormatedStructure = oldSideChainInstaller(currentCandidateFormatedBackbone, monteCarloOldStructure);
        newCandidateFormatedStructure = sideChainRotator(currentCandidateFormatedStructure, 10, sideChainRotatorResidueIndex);
        newCandidateConteactedResidueIndex = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructure, 2.351);
        
        for labelingIndex = 1:length(residue1List)
            [newCandidateSimulatedDistanceDistributionY(:, labelingIndex), ~] = DEERefineMTSSLLabeling(newCandidateFormatedStructure, residue1List(labelingIndex), residue2List(labelingIndex), R1_20210523);
        end
        [~, newCandidateDEERMaximalIndex] = max(newCandidateSimulatedDistanceDistributionY);
        monteCarloNewDistanceDistribution = newCandidateDEERMaximalIndex;
        monteCarloClashesTemperature = 3;
        [DEERIfMovingCriterion, RMSEFirstPhase(monteCarloSteps+1)] = monteCarloMetropolisCriterionGenerator(monteCarloOldDistanceDistribution, monteCarloNewDistanceDistribution, targetDEERMaximalIndex, 0.05);
        [clashesIfMovingCriterion, clashedResidueNumber(monteCarloSteps+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndex), length(newCandidateConteactedResidueIndex), 0, monteCarloClashesTemperature);
        if DEERIfMovingCriterion == 1 && clashesIfMovingCriterion == 1
            currentStatusUpdator(currentStatusFileName, 2)
            pdbSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructure, 1, monteCarloSteps); % save for test
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
                      R1LibraryMatFileName R1_20210523
        
    
    phiPsiAngleVariationAtLoopPhaseTwo = 0.1;
    clashesCriterionPhaseTwo = 2.3;
    outputTrajectoryPDBFileName = strcat(timeStamp, '_Trajectory_', num2str(1), '_', num2str(structureIndex), '.pdb');
    
    monteCarloOldContactedResidueIndexPhaseTwo = formatedPDB2contactedResidueNumbers(monteCarloOldStructurePhaseTwo, clashesCriterionPhaseTwo);
    clashedResidueNumberPhaseTwo(1) = length(monteCarloOldContactedResidueIndexPhaseTwo);
    for phaseTwoSteps = 1:30000
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
        currentCandidateFormatedStructurePhaseTwo = oldSideChainInstaller(currentCandidateFormatedBackbonePhaseTwo, monteCarloOldStructurePhaseTwo);
        newCandidateFormatedStructurePhaseTwo = sideChainRotator(currentCandidateFormatedStructurePhaseTwo, 25, monteCarloOldContactedResidueIndexPhaseTwo);
        newCandidateConteactedResidueIndexPhaseTwo = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructurePhaseTwo, clashesCriterionPhaseTwo);
        [clashesIfMovingCriterionPhaseTwo, clashedResidueNumberPhaseTwo(phaseTwoSteps+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndexPhaseTwo), length(newCandidateConteactedResidueIndexPhaseTwo), 0, 0.6);
        if clashesIfMovingCriterionPhaseTwo==1
            currentStatusUpdator(currentStatusFileName, 4)
            monteCarloOldContactedResidueIndexPhaseTwo = newCandidateConteactedResidueIndexPhaseTwo;
            monteCarloOldStructurePhaseTwo = newCandidateFormatedStructurePhaseTwo;
            monteCarloOldGeometryPhaseTwo = monteCarloNewCandidateGeometryPhaseTwo;
%             pdbSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructurePhaseTwo, 1, phaseTwoSteps); % save for test
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
                      R1LibraryMatFileName R1_20210523
    phiPsiAngleVariationAtLoopPhaseThree = 0.1;
    outputTrajectoryPDBFileName = strcat(timeStamp, '_Trajectory_', num2str(1), '_', num2str(structureIndex), '.pdb');
    
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
    monteCarloOldContactedResidueIndexPhaseThree = formatedPDB2contactedResidueNumbers(monteCarloOldStructurePhaseThree, 2.351);
    clashedResidueNumberPhaseThree(1) = length(monteCarloOldContactedResidueIndexPhaseThree);
    for monteCarloStepsPhaseThree = 1:30000
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
        currentCandidateMinimalBackboneNonBondedDistancePhaseThree = formatedBackbone2MinimalBackboneNonBondedDistance(currentCandidateFormatedBackbonePhaseThree);
        if minimalBackboneNonBondedDistance*0.98 > currentCandidateMinimalBackboneNonBondedDistancePhaseThree 
            currentStatusUpdator(currentStatusFileName, 6)
            continue
        end
        currentCandidateFormatedStructurePhaseThree = oldSideChainInstaller(currentCandidateFormatedBackbonePhaseThree, monteCarloOldStructurePhaseThree);
        newCandidateFormatedStructurePhaseThree = sideChainRotator(currentCandidateFormatedStructurePhaseThree, 5, sideChainRotatorResidueIndex);
        newCandidateConteactedResidueIndexPhaseThree = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructurePhaseThree, 2.351);
        
        for labelingIndex = 1:length(residue1List)
            [newCandidateSimulatedDistanceDistributionYPhaseThree(:, labelingIndex), ~] = DEERefineMTSSLLabeling(newCandidateFormatedStructurePhaseThree, residue1List(labelingIndex), residue2List(labelingIndex), R1_20210523);
        end
        
        [~, newCandidateDEERMaximalIndexPhaseThree] = max(newCandidateSimulatedDistanceDistributionYPhaseThree);
        monteCarloNewDistanceDistributionPhaseThree = newCandidateDEERMaximalIndexPhaseThree;
        monteCarloClashesTemperaturePhaseThree = 1;
        [DEERIfMovingCriterion, RMSEThirdPhase(monteCarloStepsPhaseThree+1)] = monteCarloMetropolisCriterionGenerator(monteCarloOldDistanceDistributionPhaseThree, monteCarloNewDistanceDistributionPhaseThree, targetDEERMaximalIndex, 0.05);
        [clashesIfMovingCriterion, clashedResidueNumberThirdPhase(monteCarloStepsPhaseThree+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndexPhaseThree), length(newCandidateConteactedResidueIndexPhaseThree), 0, monteCarloClashesTemperaturePhaseThree);
        
        if DEERIfMovingCriterion == 1 && clashesIfMovingCriterion == 1
            currentStatusUpdator(currentStatusFileName, 7)
%             pdbSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructurePhaseThree, 1, monteCarloStepsPhaseThree); % save for test
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
                      R1LibraryMatFileName R1_20210523
         
    
    
    
    
    phiPsiAngleVariationAtLoopPhaseFour = 0.1;
    clashesCriterionPhaseFour = 2.32;
    outputTrajectoryPDBFileName = strcat(timeStamp, '_Trajectory_', num2str(1), '_', num2str(structureIndex), '.pdb');
    
    monteCarloOldContactedResidueIndexPhaseFour = formatedPDB2contactedResidueNumbers(monteCarloOldStructurePhaseFour, clashesCriterionPhaseFour);
    clashedResidueNumberPhaseFour(1) = length(monteCarloOldContactedResidueIndexPhaseFour);
    for phaseFourSteps = 1:30000
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
%             pdbSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructurePhaseFour, 1, phaseFourSteps); % save for test
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
                      R1LibraryMatFileName R1_20210523
    phiPsiAngleVariationAtLoopPhaseFive = 0.1;
    outputTrajectoryPDBFileName = strcat(timeStamp, '_Trajectory_', num2str(1), '_', num2str(structureIndex), '.pdb');
    
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
    monteCarloOldContactedResidueIndexPhaseFive = formatedPDB2contactedResidueNumbers(monteCarloOldStructurePhaseFive, 2.351);
    clashedResidueNumberPhaseFive(1) = length(monteCarloOldContactedResidueIndexPhaseFive);
    
    for monteCarloStepsPhaseFive = 1:30000
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
        currentCandidateMinimalBackboneNonBondedDistancePhaseFive = formatedBackbone2MinimalBackboneNonBondedDistance(currentCandidateFormatedBackbonePhaseFive);
        if minimalBackboneNonBondedDistance*0.98 > currentCandidateMinimalBackboneNonBondedDistancePhaseFive 
            currentStatusUpdator(currentStatusFileName, 11)
            continue
        end
        currentCandidateFormatedStructurePhaseFive = oldSideChainInstaller(currentCandidateFormatedBackbonePhaseFive, monteCarloOldStructurePhaseFive);
        newCandidateFormatedStructurePhaseFive = sideChainRotator(currentCandidateFormatedStructurePhaseFive, 5, sideChainRotatorResidueIndex);
        newCandidateConteactedResidueIndexPhaseFive = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructurePhaseFive, 2.351);
        
        for labelingIndex = 1:length(residue1List)
            [newCandidateSimulatedDistanceDistributionYPhaseFive(:, labelingIndex), ~] = DEERefineMTSSLLabeling(newCandidateFormatedStructurePhaseFive, residue1List(labelingIndex), residue2List(labelingIndex), R1_20210523);
        end
        
        [~, newCandidateDEERMaximalIndexPhaseFive] = max(newCandidateSimulatedDistanceDistributionYPhaseFive);
        monteCarloNewDistanceDistributionPhaseFive = newCandidateDEERMaximalIndexPhaseFive;
        monteCarloClashesTemperaturePhaseFive = 1;
        [DEERIfMovingCriterion, RMSEFifthPhase(monteCarloStepsPhaseFive+1)] = monteCarloMetropolisCriterionGenerator(monteCarloOldDistanceDistributionPhaseFive, monteCarloNewDistanceDistributionPhaseFive, targetDEERMaximalIndex, 0.05);
        [clashesIfMovingCriterion, clashedResidueNumberFifthPhase(monteCarloStepsPhaseFive+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndexPhaseFive), length(newCandidateConteactedResidueIndexPhaseFive), 0, monteCarloClashesTemperaturePhaseFive);
        
        if DEERIfMovingCriterion == 1 && clashesIfMovingCriterion == 1
            currentStatusUpdator(currentStatusFileName, 12)
%             pdbSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructurePhaseFive, 1, monteCarloStepsPhaseFive); % save for test
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
                      R1LibraryMatFileName R1_20210523
    
    
    
    
    phiPsiAngleVariationAtLoopPhaseSix = 0.1;
    clashesCriterionPhaseSix = 2.34;
    outputTrajectoryPDBFileName = strcat(timeStamp, '_Trajectory_', num2str(1), '_', num2str(structureIndex), '.pdb');
    
    monteCarloOldContactedResidueIndexPhaseSix = formatedPDB2contactedResidueNumbers(monteCarloOldStructurePhaseSix, clashesCriterionPhaseSix);
    clashedResidueNumberPhaseSix(1) = length(monteCarloOldContactedResidueIndexPhaseSix);
    for phaseSixSteps = 1:30000
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
        currentCandidateFormatedStructurePhaseSix = oldSideChainInstaller(currentCandidateFormatedBackbonePhaseSix, monteCarloOldStructurePhaseSix);
        newCandidateFormatedStructurePhaseSix = sideChainRotator(currentCandidateFormatedStructurePhaseSix, 25, monteCarloOldContactedResidueIndexPhaseSix);
        newCandidateConteactedResidueIndexPhaseSix = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructurePhaseSix, clashesCriterionPhaseSix);
        [clashesIfMovingCriterionPhaseSix, clashedResidueNumberPhaseSix(phaseSixSteps+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndexPhaseSix), length(newCandidateConteactedResidueIndexPhaseSix), 0, 0.6);
        if clashesIfMovingCriterionPhaseSix==1
            currentStatusUpdator(currentStatusFileName, 14)
            monteCarloOldContactedResidueIndexPhaseSix = newCandidateConteactedResidueIndexPhaseSix;
            monteCarloOldStructurePhaseSix = newCandidateFormatedStructurePhaseSix;
            monteCarloOldGeometryPhaseSix = monteCarloNewCandidateGeometryPhaseSix;
%             pdbSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructurePhaseSix, 1, phaseSixSteps); % save for test
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
                      R1LibraryMatFileName R1_20210523
    phiPsiAngleVariationAtLoopPhaseSeven = 0.1;
    outputTrajectoryPDBFileName = strcat(timeStamp, '_Trajectory_', num2str(1), '_', num2str(structureIndex), '.pdb');
    
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
    monteCarloOldContactedResidueIndexPhaseSeven = formatedPDB2contactedResidueNumbers(monteCarloOldStructurePhaseSeven, 2.351);
    clashedResidueNumberPhaseSeven(1) = length(monteCarloOldContactedResidueIndexPhaseSeven);
    
    for monteCarloStepsPhaseSeven = 1:30000
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
        currentCandidateMinimalBackboneNonBondedDistancePhaseSeven = formatedBackbone2MinimalBackboneNonBondedDistance(currentCandidateFormatedBackbonePhaseSeven);
        if minimalBackboneNonBondedDistance*0.98 > currentCandidateMinimalBackboneNonBondedDistancePhaseSeven 
            currentStatusUpdator(currentStatusFileName, 16)
            continue
        end
        currentCandidateFormatedStructurePhaseSeven = oldSideChainInstaller(currentCandidateFormatedBackbonePhaseSeven, monteCarloOldStructurePhaseSeven);
        newCandidateFormatedStructurePhaseSeven = sideChainRotator(currentCandidateFormatedStructurePhaseSeven, 5, sideChainRotatorResidueIndex);
        newCandidateConteactedResidueIndexPhaseSeven = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructurePhaseSeven, 2.351);
        
        for labelingIndex = 1:length(residue1List)
            [newCandidateSimulatedDistanceDistributionYPhaseSeven(:, labelingIndex), ~] = DEERefineMTSSLLabeling(newCandidateFormatedStructurePhaseSeven, residue1List(labelingIndex), residue2List(labelingIndex), R1_20210523);
        end
        
        [~, newCandidateDEERMaximalIndexPhaseSeven] = max(newCandidateSimulatedDistanceDistributionYPhaseSeven);
        monteCarloNewDistanceDistributionPhaseSeven = newCandidateDEERMaximalIndexPhaseSeven;
        monteCarloClashesTemperaturePhaseSeven = 1;
        [DEERIfMovingCriterion, RMSESeventhPhase(monteCarloStepsPhaseSeven+1)] = monteCarloMetropolisCriterionGenerator(monteCarloOldDistanceDistributionPhaseSeven, monteCarloNewDistanceDistributionPhaseSeven, targetDEERMaximalIndex, 0.05);
        [clashesIfMovingCriterion, clashedResidueNumberSeventhPhase(monteCarloStepsPhaseSeven+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndexPhaseSeven), length(newCandidateConteactedResidueIndexPhaseSeven), 0, monteCarloClashesTemperaturePhaseSeven);
        
        if DEERIfMovingCriterion == 1 && clashesIfMovingCriterion == 1
            currentStatusUpdator(currentStatusFileName, 17)
%             pdbSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructurePhaseSeven, 1, monteCarloStepsPhaseSeven); % save for test
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
                      R1LibraryMatFileName R1_20210523
    
    
    
    
    phiPsiAngleVariationAtLoopPhaseEight = 0.1;
    clashesCriterionPhaseEight = 2.35;
    outputTrajectoryPDBFileName = strcat(timeStamp, '_Trajectory_', num2str(1), '_', num2str(structureIndex), '.pdb');
    
    monteCarloOldContactedResidueIndexPhaseEight = formatedPDB2contactedResidueNumbers(monteCarloOldStructurePhaseEight, clashesCriterionPhaseEight);
    clashedResidueNumberPhaseEight(1) = length(monteCarloOldContactedResidueIndexPhaseEight);
    for phaseEightSteps = 1:30000
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
        currentCandidateFormatedStructurePhaseEight = oldSideChainInstaller(currentCandidateFormatedBackbonePhaseEight, monteCarloOldStructurePhaseEight);
        newCandidateFormatedStructurePhaseEight = sideChainRotator(currentCandidateFormatedStructurePhaseEight, 25, monteCarloOldContactedResidueIndexPhaseEight);
        newCandidateConteactedResidueIndexPhaseEight = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructurePhaseEight, clashesCriterionPhaseEight);
        [clashesIfMovingCriterionPhaseEight, clashedResidueNumberPhaseEight(phaseEightSteps+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndexPhaseEight), length(newCandidateConteactedResidueIndexPhaseEight), 0, 0.6);
        if clashesIfMovingCriterionPhaseEight==1
            currentStatusUpdator(currentStatusFileName, 19)
            monteCarloOldContactedResidueIndexPhaseEight = newCandidateConteactedResidueIndexPhaseEight;
            monteCarloOldStructurePhaseEight = newCandidateFormatedStructurePhaseEight;
            monteCarloOldGeometryPhaseEight = monteCarloNewCandidateGeometryPhaseEight;
%             pdbSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructurePhaseEight, 1, phaseEightSteps); % save for test
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
                      R1LibraryMatFileName R1_20210523
    phiPsiAngleVariationAtLoopPhaseNine = 0.1;
    outputTrajectoryPDBFileName = strcat(timeStamp, '_Trajectory_', num2str(1), '_', num2str(structureIndex), '.pdb');
    
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
    monteCarloOldContactedResidueIndexPhaseNine = formatedPDB2contactedResidueNumbers(monteCarloOldStructurePhaseNine, 2.351);
    clashedResidueNumberPhaseNine(1) = length(monteCarloOldContactedResidueIndexPhaseNine);
    
    for monteCarloStepsPhaseNine = 1:30000
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
        currentCandidateMinimalBackboneNonBondedDistancePhaseNine = formatedBackbone2MinimalBackboneNonBondedDistance(currentCandidateFormatedBackbonePhaseNine);
        if minimalBackboneNonBondedDistance*0.98 > currentCandidateMinimalBackboneNonBondedDistancePhaseNine 
            currentStatusUpdator(currentStatusFileName, 21)
            continue
        end
        currentCandidateFormatedStructurePhaseNine = oldSideChainInstaller(currentCandidateFormatedBackbonePhaseNine, monteCarloOldStructurePhaseNine);
        newCandidateFormatedStructurePhaseNine = sideChainRotator(currentCandidateFormatedStructurePhaseNine, 5, sideChainRotatorResidueIndex);
        newCandidateConteactedResidueIndexPhaseNine = formatedPDB2contactedResidueNumbers(newCandidateFormatedStructurePhaseNine, 2.351);
        
        for labelingIndex = 1:length(residue1List)
            [newCandidateSimulatedDistanceDistributionYPhaseNine(:, labelingIndex), ~] = DEERefineMTSSLLabeling(newCandidateFormatedStructurePhaseNine, residue1List(labelingIndex), residue2List(labelingIndex), R1_20210523);
        end
        
        [~, newCandidateDEERMaximalIndexPhaseNine] = max(newCandidateSimulatedDistanceDistributionYPhaseNine);
        monteCarloNewDistanceDistributionPhaseNine = newCandidateDEERMaximalIndexPhaseNine;
        monteCarloClashesTemperaturePhaseNine = 1;
        [DEERIfMovingCriterion, RMSENinethPhase(monteCarloStepsPhaseNine+1)] = monteCarloMetropolisCriterionGenerator(monteCarloOldDistanceDistributionPhaseNine, monteCarloNewDistanceDistributionPhaseNine, targetDEERMaximalIndex, 0.05);
        [clashesIfMovingCriterion, clashedResidueNumberNinethPhase(monteCarloStepsPhaseNine+1)] = monteCarloMetropolisCriterionGenerator(length(monteCarloOldContactedResidueIndexPhaseNine), length(newCandidateConteactedResidueIndexPhaseNine), 0, monteCarloClashesTemperaturePhaseNine);
        
        if DEERIfMovingCriterion == 1 && clashesIfMovingCriterion == 1
            currentStatusUpdator(currentStatusFileName, 22)
%             pdbSaver(outputTrajectoryPDBFileName, newCandidateFormatedStructurePhaseNine, 1, monteCarloStepsPhaseNine); % save for test
            monteCarloOldContactedResidueIndexPhaseNine = newCandidateConteactedResidueIndexPhaseNine;
            monteCarloOldStructurePhaseNine = newCandidateFormatedStructurePhaseNine;
            monteCarloOldDistanceDistributionPhaseNine = monteCarloNewDistanceDistributionPhaseNine;
            monteCarloOldGeometryPhaseNine = monteCarloNewCandidateGeometryPhaseNine;
            if RMSENinethPhase(monteCarloStepsPhaseNine+1) <= RMSE+0.2 && clashedResidueNumberNinethPhase(monteCarloStepsPhaseNine+1) <= maximalClashes+2
                increaseFifthRMSEPassedNumber(fifthRMSEPassedNumberFileName)

                fifthRMSEPassedNumber = load(fifthRMSEPassedNumberFileName);
                candidateGenerator(fifthRMSEPassedNumber, newCandidateFormatedStructurePhaseNine, timeStamp)

                currentStatusUpdator(currentStatusFileName, 23)

                if fifthRMSEPassedNumber == numberOfStructure
                    FINALPDBENSEMBLEFileGenerator(timeStamp)
                    pause(5)
                    FINALPDBENSEMBLEFile = dir('*FINALPDBENSEMBLE_*');
                    FINALPDBENSEMBLEFileName = FINALPDBENSEMBLEFile.name;

    %                 FINALPDBENSEMBLE = pdbModelsLoader(FINALPDBENSEMBLEFileName);
    %                 for FINALPDBNumber = 1:length(FINALPDBENSEMBLE(1, 1, :))
    %                     for labelingIndex = 1:length(residue1List)
    %                         [distanceDistributionFinal(:, labelingIndex, FINALPDBNumber), ~] = DEERefineMTSSLLabeling(FINALPDBENSEMBLE(:, :, FINALPDBNumber), residue1List(labelingIndex), residue2List(labelingIndex));
    %                         distanceDistributionFinal(:, labelingIndex, FINALPDBNumber) = distanceDistributionFinal(:, labelingIndex, FINALPDBNumber)/sum(distanceDistributionFinal(:, labelingIndex, FINALPDBNumber));
    %                     end
    %                 end
    %                 for FINALPDBNumber = 1:length(FINALPDBENSEMBLE(1, 1, :))
    %                     distanceDistributionFinalForJSD(:, FINALPDBNumber) = reshape(distanceDistributionFinal(:, :, FINALPDBNumber), [], 1);
    %                 end
                    [distanceDistributionFinalForJSD, FINALPDBENSEMBLE] = distanceDistributionFinalForJSDGenerator(FINALPDBENSEMBLEFileName, residue1List, residue2List);
    %                 for distanceDistributionNumber = 1:length(distanceDistributionFullPath)
    %                     targetDistanceDistributionFinal(:, :, distanceDistributionNumber) = distanceDistributionFinalLoader(distanceDistributionFullPath(distanceDistributionNumber));
    %                 end
    %                 targetDistanceDistributionFinalForJSD = reshape(targetDistanceDistributionFinal(:, 2, :), [], 1);
                    [targetDistanceDistributionFinalForJSD, targetDistanceDistributionFinal] = targetDistanceDistributionFinalForJSDGenerator(distanceDistributionFullPath);



                    for FINALPDBNumber = 1:length(FINALPDBENSEMBLE(1, 1, :))
                        JensenShannonDivergence(FINALPDBNumber) = convertingDistanceDistributions2JensenShannonDivergence(targetDistanceDistributionFinalForJSD, distanceDistributionFinalForJSD(:, FINALPDBNumber));
                    end

    %                 FINALPDBFileName = strcat("FINALPDB_",  timeStamp, ".pdb");
    %                 [~, minimalJSDIndex] = min(JensenShannonDivergence);
    %                 pdbSaver(FINALPDBFileName, FINALPDBENSEMBLE(:, :, minimalJSDIndex));
                    minimalJSDIndex = minimalJSDPDBGenerator(FINALPDBENSEMBLE, JensenShannonDivergence, timeStamp);
    %                 for labelingIndex = 1:length(residue1List)
    %                     [DistanceDistributionFinal(:, 2, labelingIndex), DistanceDistributionFinal(:, 1, labelingIndex)] = DEERefineMTSSLLabeling(FINALPDBENSEMBLE(:, :, minimalJSDIndex), residue1List(labelingIndex), residue2List(labelingIndex));
    %                 end
    %                 
    %                 normalizedTargetDistanceDistributionFinal = targetDistanceDistributionFinal;
    %                 for distanceDistributionNumber = 1:length(distanceDistributionFullPath)
    %                     normalizedTargetDistanceDistributionFinal(:, 2, distanceDistributionNumber) = targetDistanceDistributionFinal(:, 2, distanceDistributionNumber)/max(targetDistanceDistributionFinal(:, 2, distanceDistributionNumber));
    %                 end
    %                 
    %                 DEERefineFinalDistanceDistributionFileName = strcat("FINALPr_",  timeStamp, ".mat");
    %                 DEERefineFinalDistanceDistribution = struct('targetDistanceDistributionFinal', normalizedTargetDistanceDistributionFinal,...
    %                                                             'DistanceDistributionFinal', DistanceDistributionFinal);
    %                 save(DEERefineFinalDistanceDistributionFileName, '-struct', 'DEERefineFinalDistanceDistribution');
                    DEERefineFinalFileGenerator(residue1List, residue2List, FINALPDBENSEMBLE, minimalJSDIndex, targetDistanceDistributionFinal, distanceDistributionFullPath, timeStamp)
    %                 copyfile(DEERefineFinalDistanceDistributionFileName, '..');
    %                 copyfile(FINALPDBENSEMBLEFileName, '..');
    %                 copyfile(FINALPDBFileName, '..');
                    DEERefineFinalFileCopier(timeStamp)
                    closereq
                    return
                elseif fifthRMSEPassedNumber > numberOfStructure
                    return
                end
                % GLOBAL X
                % X = 0;
                % X = X+1;
                % finalSturcture(:, :, X) =
                % newCandidateFormatedStructurePhaseNine;
                % if X >numberOfStructure
                %   cancel(job)
                % end
    %             pdbSaver(phaseOnePassedFileName, newCandidateFormatedStructure);
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
