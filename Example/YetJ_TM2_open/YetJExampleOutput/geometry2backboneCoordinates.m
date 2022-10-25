function backboneCoordinates = geometry2backboneCoordinates(geometry)  
    geometry.length.backboneCoordinatesTmp = reshape([geometry.length.NCA  geometry.length.CAC [geometry.length.CN1; 0]]', 1, length(geometry.length.NCA)*3)';
    geometry.length.backboneCoordinatesTmp(end) = [];
    geometry.length.backbone = geometry.length.backboneCoordinatesTmp;
    geometry.length = rmfield(geometry.length, 'backboneCoordinatesTmp');

    geometry.angle.backboneCoordinatesTmp = reshape([geometry.angle.NCAC  [geometry.angle.CACN1; 0] [geometry.angle.CN1CA1; 0]]', 1, length(geometry.angle.NCAC)*3)';
    geometry.angle.backboneCoordinatesTmp(end-1:end) = [];
    geometry.angle.backbone = geometry.angle.backboneCoordinatesTmp;
    geometry.angle = rmfield(geometry.angle, 'backboneCoordinatesTmp');

    geometry.dihedral.backboneCoordinatesTmp = reshape([[geometry.dihedral.NCACN1]  [geometry.dihedral.CACN1CA1] [geometry.dihedral.CN1CA1C1]]', 1, length(geometry.dihedral.NCACN1)*3)';
    geometry.dihedral.backbone = geometry.dihedral.backboneCoordinatesTmp;
    geometry.dihedral = rmfield(geometry.dihedral, 'backboneCoordinatesTmp');
    backboneCoordinatesTmp = backboneGeometry2backboneCoordinates(geometry.length.backbone, geometry.angle.backbone, geometry.dihedral.backbone);
    for i = 1:length(geometry.length.NCA)
        bondLength(:, i) = [geometry.length.NCA(i); geometry.length.CAC(i); geometry.length.CO(i)];
        bondAngle(:, i) = [geometry.angle.NCAC(i); geometry.angle.CACO(i)];
        bondDihedral(:, i) = geometry.dihedral.NCACO(i);
        backboneCoordinatesBeforeRotations(:, :, i) = backboneGeometry2backboneCoordinates(bondLength(:, i), bondAngle(:, i), bondDihedral(:, i));
        [firstRotatingMatrix(:, :, i), secondRotatingMatrix(:, :, i), thirdRotatingMatrix(:, :, i)] = rotatingMatricesGenerator(backboneCoordinatesTmp, i);
        backboneCoordinatesAfterRotations(:, :, i) = backboneCoordinatesBeforeRotations(:, :, i)*transpose(thirdRotatingMatrix(:, :, i))*transpose(secondRotatingMatrix(:, :, i))*transpose(firstRotatingMatrix(:, :, i))+backboneCoordinatesTmp(3*i-2, :);
    end
    backboneCoordinates = reshape(pagetranspose(backboneCoordinatesAfterRotations), 3, 4*length(geometry.length.NCA))';
    
    function backboneCoordinates = backboneGeometry2backboneCoordinates(backboneLength, backboneAngle, backboneDihedral)
        flippedBackboneLength = flip(backboneLength);
        extendedBackboneAngle = 180-backboneAngle;
        flippedExtendedBackboneAngle = flip(extendedBackboneAngle);
        extendedBackboneDihedral = [0; backboneDihedral];
        flippedExtendedBackboneDihedral = flip(extendedBackboneDihedral);
        backboneCoordinates = [0 0 0];
        for k = 1:length(flippedExtendedBackboneDihedral)
            backboneCoordinates = backboneCoordinates+[0 0 flippedBackboneLength(k)];
            backboneCoordinates = backboneCoordinates*rotatingMatrixGenerator([1, 0, 0], flippedExtendedBackboneAngle(k)/180*pi);
            backboneCoordinates = backboneCoordinates*rotatingMatrixGenerator([0, 0, 1], flippedExtendedBackboneDihedral(k)/180*pi);
            backboneCoordinates = [0 0 0; backboneCoordinates];
        end
        backboneCoordinates = backboneCoordinates+[0 0 flippedBackboneLength(end)];
        backboneCoordinates = [0 0 0; backboneCoordinates];
    end
    function [firstRotatingMatrix, secondRotatingMatrix, thirdRotatingMatrix] = rotatingMatricesGenerator(backboneCoordinates, residueIndex)
        originBackboneCoordinates = backboneCoordinates-backboneCoordinates(residueIndex*3-2, :);
        firstRotatingMatrix = rotatingMatrixGenerator([1 1 1], 0.001);
        randomlyRotatedOriginBackboneCoordinates = originBackboneCoordinates*firstRotatingMatrix;
        secondRotatingMatrix = rotatingMatrixGenerator(cross(randomlyRotatedOriginBackboneCoordinates(residueIndex*3-1, :), [0 0 1]), acos(dot(randomlyRotatedOriginBackboneCoordinates(residueIndex*3-1, :)/norm(randomlyRotatedOriginBackboneCoordinates(residueIndex*3-1, :)), [0 0 1], 2)));
        alignedBackboneCoordinates = randomlyRotatedOriginBackboneCoordinates*secondRotatingMatrix;
        thirdRotatingMatrix = rotatingMatrixGenerator([0 0 1], sign(sum(cross([alignedBackboneCoordinates(residueIndex*3, 1:2) 0], [0 -1 0])))*acos(dot([alignedBackboneCoordinates(residueIndex*3, 1:2) 0]/norm([alignedBackboneCoordinates(residueIndex*3, 1:2) 0]), [0 -1 0], 2)));
    end
    
end