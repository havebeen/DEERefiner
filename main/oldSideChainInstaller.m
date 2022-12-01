function oldSideChainInstalledFormatedStructure = oldSideChainInstaller(newFormatedBackbone, oldFormatedStructure)
    residueIndex = unique(double(oldFormatedStructure(:, 7)));
    oldSideChainCoordinates = nan(14, 3, length(residueIndex));
    oldFormatedSideChain = strings(14, 16, length(residueIndex));
    newBackboneCoordinates = nan(4, 3, length(residueIndex));
    pagedNewFormatedBackbone = strings(4, 16, length(residueIndex));
    
    residueIndexString = strrep(oldFormatedStructure(:, 7), ' ', '');
    backboneResidueIndexString = strrep(oldFormatedStructure(oldFormatedStructure(:, 3)==" N  "|...
                                   oldFormatedStructure(:, 3)==" CA "|...
                                   oldFormatedStructure(:, 3)==" C  "|...
                                   oldFormatedStructure(:, 3)==" O  ", 7), ' ', '');
    
    for loopResidueIndex = 1:length(residueIndex)
    oldSideChainCoordinates(1:length(oldFormatedStructure(residueIndexString==string(residueIndex(loopResidueIndex)), 1)), :, loopResidueIndex) = ...
        double(oldFormatedStructure(residueIndexString == string(residueIndex(loopResidueIndex)), 9:11));
    oldFormatedSideChain(1:length(oldFormatedStructure(residueIndexString==string(residueIndex(loopResidueIndex)), 1)), :, loopResidueIndex) = ...
        oldFormatedStructure(residueIndexString == string(residueIndex(loopResidueIndex)), :);
    
    newBackboneCoordinates(:, :, loopResidueIndex) = double(newFormatedBackbone(backboneResidueIndexString == string(residueIndex(loopResidueIndex)), 9:11));
    pagedNewFormatedBackbone(:, :, loopResidueIndex) = newFormatedBackbone(backboneResidueIndexString == string(residueIndex(loopResidueIndex)), :);
    end
    
    newFormatedBackboneTargetCACoordinates = pagetranspose(reshape(double(newFormatedBackbone(newFormatedBackbone(:, 3) == " CA ", 9:11))', 3, 1, []));
    newFormatedBackboneProteinCoordinatesFirstMoved = newBackboneCoordinates-newFormatedBackboneTargetCACoordinates;
    newFormatedBackboneTargetCCoordinatesFirstMoved = pagetranspose(reshape(double(newFormatedBackbone(newFormatedBackbone(:, 3) == " C  ", 9:11))', 3, 1, []))-newFormatedBackboneTargetCACoordinates;
    newFormatedBackboneVectorCACFirstMoved = newFormatedBackboneTargetCCoordinatesFirstMoved;
    newFormatedBackboneCrossCACzFirstMoved = cross(newFormatedBackboneVectorCACFirstMoved, repmat([0 0 1], 1, 1, length(residueIndex)));
    normalizedNewFormatedBackboneVectorCACFirstMoved = sqrt(sum(newFormatedBackboneVectorCACFirstMoved.^2));
    newFormatedBackboneAngleCACzFirstMoved = acos(dot(newFormatedBackboneVectorCACFirstMoved, repmat([0 0 1], 1, 1, length(residueIndex)))./normalizedNewFormatedBackboneVectorCACFirstMoved/1);
    
    newFormatedBackboneTmp = newFormatedBackbone;
    newFormatedBackboneProteinCoordinatesSecondMoved = pagemtimes(newFormatedBackboneProteinCoordinatesFirstMoved,rotatingMatrixGenerator(newFormatedBackboneCrossCACzFirstMoved, newFormatedBackboneAngleCACzFirstMoved));
    newFormatedBackboneTmp(:, 9:11) = string(reshape(pagetranspose(newFormatedBackboneProteinCoordinatesSecondMoved), 3, [])');
    newFormatedBackboneTargetNCoordinatesSecondMoved = reshape(double(newFormatedBackboneTmp(newFormatedBackboneTmp(:, 3) == " N  ", 9:11))', 1, 3, []);
    newFormatedBackboneTargetCCoordinatesSecondMoved = reshape(double(newFormatedBackboneTmp(newFormatedBackboneTmp(:, 3) == " C  ", 9:11))', 1, 3, []);
    newFormatedBackboneVectorCANSecondMoved = newFormatedBackboneTargetNCoordinatesSecondMoved;
    projectedNewFormatedBackboneVectorCANSecondMoved = [newFormatedBackboneVectorCANSecondMoved(:, 1:2, :) zeros(1, 1, length(residueIndex))];
    newFormatedBackboneSecondMovedRotatingSign = sign(dot(cross(projectedNewFormatedBackboneVectorCANSecondMoved, repmat([1 0 0], 1, 1, length(residueIndex))), newFormatedBackboneTargetCCoordinatesSecondMoved));
    normalizedProjectedNewFormatedBackboneVectorCANSecondMoved = sqrt(sum(projectedNewFormatedBackboneVectorCANSecondMoved.^2));
    newFormatedBackboneAngleCANxSecondMoved = acos(dot(projectedNewFormatedBackboneVectorCANSecondMoved, repmat([1 0 0], 1, 1, length(residueIndex)))./normalizedProjectedNewFormatedBackboneVectorCANSecondMoved/1);
    newFormatedBackboneProteinCoordinatesThirdMoved = pagemtimes(newFormatedBackboneProteinCoordinatesSecondMoved, rotatingMatrixGenerator(newFormatedBackboneTargetCCoordinatesSecondMoved, newFormatedBackboneSecondMovedRotatingSign.*newFormatedBackboneAngleCANxSecondMoved));

    
    
    oldFormatedStructureTargetCACoordinates = pagetranspose(reshape(double(oldFormatedStructure(oldFormatedStructure(:, 3) == " CA ", 9:11))', 3, 1, []));
    oldFormatedStructureProteinCoordinatesFirstMoved = oldSideChainCoordinates-oldFormatedStructureTargetCACoordinates;
    oldFormatedStructureTargetCCoordinatesFirstMoved = pagetranspose(reshape(double(oldFormatedStructure(oldFormatedStructure(:, 3) == " C  ", 9:11))', 3, 1, []))-oldFormatedStructureTargetCACoordinates;
    oldFormatedStructureVectorCACFirstMoved = oldFormatedStructureTargetCCoordinatesFirstMoved;
    oldFormatedStructureCrossCACzFirstMoved = cross(oldFormatedStructureVectorCACFirstMoved, repmat([0 0 1], 1, 1, length(residueIndex)));
    normalizedOldFormatedStructureVectorCACFirstMoved = sqrt(sum(oldFormatedStructureVectorCACFirstMoved.^2));
    oldFormatedStructureAngleCACzFirstMoved = acos(dot(oldFormatedStructureVectorCACFirstMoved, repmat([0 0 1], 1, 1, length(residueIndex)))./normalizedOldFormatedStructureVectorCACFirstMoved/1);
    
    
    oldFormatedStructureTmp = oldFormatedStructure;
    oldFormatedStructureProteinCoordinatesSecondMoved = pagemtimes(oldFormatedStructureProteinCoordinatesFirstMoved,rotatingMatrixGenerator(oldFormatedStructureCrossCACzFirstMoved, oldFormatedStructureAngleCACzFirstMoved));
    reshapedOldFormatedStructureProteinCoordinatesSecondMoved = reshape(pagetranspose(oldFormatedStructureProteinCoordinatesSecondMoved), 3, [])';
    nanRemovedOldFormatedStructureProteinCoordinatesSecondMoved = rmmissing(reshapedOldFormatedStructureProteinCoordinatesSecondMoved);
    oldFormatedStructureTmp(:, 9:11) = string(nanRemovedOldFormatedStructureProteinCoordinatesSecondMoved);
    oldFormatedStructureTargetNCoordinatesSecondMoved = reshape(double(oldFormatedStructureTmp(oldFormatedStructureTmp(:, 3) == " N  ", 9:11))', 1, 3, []);
    oldFormatedStructureTargetCCoordinatesSecondMoved = reshape(double(oldFormatedStructureTmp(oldFormatedStructureTmp(:, 3) == " C  ", 9:11))', 1, 3, []);
    oldFormatedStructureVectorCANSecondMoved = oldFormatedStructureTargetNCoordinatesSecondMoved;
    projectedOldFormatedStructureVectorCANSecondMoved = [oldFormatedStructureVectorCANSecondMoved(:, 1:2, :) zeros(1, 1, length(residueIndex))];
    oldFormatedStructureSecondMovedRotatingSign = sign(dot(cross(projectedOldFormatedStructureVectorCANSecondMoved, repmat([1 0 0], 1, 1, length(residueIndex))), oldFormatedStructureTargetCCoordinatesSecondMoved));
    normalizedProjectedOldFormatedStructureVectorCANSecondMoved = sqrt(sum(projectedOldFormatedStructureVectorCANSecondMoved.^2));
    oldFormatedStructureAngleCANxSecondMoved = acos(dot(projectedOldFormatedStructureVectorCANSecondMoved, repmat([1 0 0], 1, 1, length(residueIndex)))./normalizedProjectedOldFormatedStructureVectorCANSecondMoved/1);
    oldFormatedStructureProteinCoordinatesThirdMoved = pagemtimes(oldFormatedStructureProteinCoordinatesSecondMoved, rotatingMatrixGenerator(oldFormatedStructureTargetCCoordinatesSecondMoved, oldFormatedStructureSecondMovedRotatingSign.*oldFormatedStructureAngleCANxSecondMoved));
    
    
    sideChainFirstMoved = pagemtimes(oldFormatedStructureProteinCoordinatesThirdMoved,rotatingMatrixGenerator(newFormatedBackboneTargetCCoordinatesSecondMoved, -1*newFormatedBackboneSecondMovedRotatingSign.*newFormatedBackboneAngleCANxSecondMoved));
    sideChainSecondMoved = pagemtimes(sideChainFirstMoved, rotatingMatrixGenerator(newFormatedBackboneCrossCACzFirstMoved, -1*newFormatedBackboneAngleCACzFirstMoved));
    sideChainThirdMoved = sideChainSecondMoved + newFormatedBackboneTargetCACoordinates;
    reshapedSideChainThirdMoved = reshape(pagetranspose(sideChainThirdMoved), 3, [])';
    sideChainThirdMovedCoordinates = rmmissing(reshapedSideChainThirdMoved);
    oldFormatedStructureTmp(:, 9:11) = string(sideChainThirdMovedCoordinates);
    oldSideChainInstalledFormatedStructure = oldFormatedStructureTmp;
    
    
    
    
    
end