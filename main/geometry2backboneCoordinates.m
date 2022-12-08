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
    %%
    originBackboneCoordinates_aa = repmat(backboneCoordinatesTmp, 1, 1, length(geometry.length.NCA));

    CA3D = reshape(backboneCoordinatesTmp(1:3:length(backboneCoordinatesTmp(:, 1)), :)', 1, 3, []);
    originBackboneCoordinates_aa = originBackboneCoordinates_aa-CA3D;

    randomlyRotatedOriginBackboneCoordinates = pagemtimes(originBackboneCoordinates_aa, repmat([1 0.0006 -0.0006; -0.0006 1 0.0006;0.0006 -0.0006 1], 1, 1, length(geometry.length.NCA)));



    for i = 1:length(CA3D(1, 1, :))
    firstCrossVector(1, :, i) = randomlyRotatedOriginBackboneCoordinates(i*3-1, :, i);
    end
    a = cross(firstCrossVector, repmat([0 0 1], 1, 1, length(CA3D(1, 1, :))));
    b = acos(dot(firstCrossVector./DNorm2(firstCrossVector), repmat([0 0 1], 1, 1, length(CA3D(1, 1, :)))));
    secondRotatingMatrix = rotatingMatrixGenerator(a, b);


    alignedBackboneCoordinates = pagemtimes(randomlyRotatedOriginBackboneCoordinates, secondRotatingMatrix);

    for i = 1:length(CA3D(1, 1, :))
        secondCrossVector(1, :, i) = alignedBackboneCoordinates(i*3, :, i);
    end
    secondCrossVector(:, 3, :) = 0;

    cross(secondCrossVector, repmat([0 -1 0], 1, 1, length(CA3D(1, 1, :))));

    c = sign(sum(cross(secondCrossVector, repmat([0 -1 0], 1, 1, length(CA3D(1, 1, :))))));




    d = acos(dot(secondCrossVector./DNorm2(secondCrossVector), repmat([0 -1 0], 1, 1, length(CA3D(1, 1, :))), 2));

    thirdRotatingMatrix = rotatingMatrixGenerator(repmat([0 0 1], 1, 1, length(CA3D(1, 1, :))), pagemtimes(c, d));
    allRotatingMatrix = pagemtimes(pagetranspose(thirdRotatingMatrix), pagetranspose(secondRotatingMatrix));
    %%
    for i = 1:length(geometry.length.NCA)
        bondLength(:, i) = [geometry.length.NCA(i); geometry.length.CAC(i); geometry.length.CO(i)];
        bondAngle(:, i) = [geometry.angle.NCAC(i); geometry.angle.CACO(i)];
        bondDihedral(:, i) = geometry.dihedral.NCACO(i);
        backboneCoordinatesBeforeRotations(:, :, i) = backboneGeometry2backboneCoordinates(bondLength(:, i), bondAngle(:, i), bondDihedral(:, i));
%         [secondRotatingMatrix(:, :, i), thirdRotatingMatrix(:, :, i)] = rotatingMatricesGenerator(backboneCoordinatesTmp, i);
        backboneCoordinatesAfterRotations(:, :, i) = backboneCoordinatesBeforeRotations(:, :, i)*allRotatingMatrix(:, :, i)*transpose([1 0.0006 -0.0006; -0.0006 1 0.0006;0.0006 -0.0006 1])+backboneCoordinatesTmp(3*i-2, :);
    end
    backboneCoordinates = reshape(pagetranspose(backboneCoordinatesAfterRotations), 3, 4*length(geometry.length.NCA))';
    

    
end
function backboneCoordinates = backboneGeometry2backboneCoordinates(backboneLength, backboneAngle, backboneDihedral)
    flippedBackboneLength = flip(backboneLength);
    extendedBackboneAngle = 180-backboneAngle;
    flippedExtendedBackboneAngle = flip(extendedBackboneAngle);
    extendedBackboneDihedral = [0; backboneDihedral];
    flippedExtendedBackboneDihedral = flip(extendedBackboneDihedral);
    backboneCoordinates = [0 0 0];
    flippedExtendedBackboneAngle3D = reshape(flippedExtendedBackboneAngle/180*pi, 1, 1, []);
    flippedExtendedBackboneDihedral3D = reshape(flippedExtendedBackboneDihedral/180*pi, 1, 1, []);
    firstRotatingMatrixWithinFunction = rotatingMatrixGenerator(repmat([1, 0, 0], 1, 1, length(flippedExtendedBackboneAngle3D)), flippedExtendedBackboneAngle3D);
    secondRotatingMatrixWithinFunction = rotatingMatrixGenerator(repmat([0 0 1], 1, 1, length(flippedExtendedBackboneDihedral3D)), flippedExtendedBackboneDihedral3D);
    rotatingMatrixWithingFunction = pagemtimes(firstRotatingMatrixWithinFunction, secondRotatingMatrixWithinFunction);
    for k = 1:length(flippedExtendedBackboneDihedral)
        backboneCoordinates = (backboneCoordinates+[0 0 flippedBackboneLength(k)])*rotatingMatrixWithingFunction(:, :, k);
%             backboneCoordinates = backboneCoordinates*firstRotatingMatrixWithinFunction(:, :, k);
%             backboneCoordinates = backboneCoordinates*secondRotatingMatrixWithinFunction(:, :, k);
        backboneCoordinates = [0 0 0; backboneCoordinates];
    end
    backboneCoordinates = backboneCoordinates+[0 0 flippedBackboneLength(end)];
    backboneCoordinates = [0 0 0; backboneCoordinates];
end
% function [secondRotatingMatrix, thirdRotatingMatrix] = rotatingMatricesGenerator(backboneCoordinates, residueIndex)
%     originBackboneCoordinates = backboneCoordinates-backboneCoordinates(residueIndex*3-2, :);
%     randomlyRotatedOriginBackboneCoordinates = originBackboneCoordinates*[1 0.0006 -0.0006; -0.0006 1 0.0006;0.0006 -0.0006 1];
%     secondRotatingMatrix = rotatingMatrixGenerator(cross(randomlyRotatedOriginBackboneCoordinates(residueIndex*3-1, :), [0 0 1]), acos(dot(randomlyRotatedOriginBackboneCoordinates(residueIndex*3-1, :)/norm(randomlyRotatedOriginBackboneCoordinates(residueIndex*3-1, :)), [0 0 1], 2)));
%     alignedBackboneCoordinates = randomlyRotatedOriginBackboneCoordinates*secondRotatingMatrix;
%     thirdRotatingMatrix = rotatingMatrixGenerator([0 0 1], sign(sum(cross([alignedBackboneCoordinates(residueIndex*3, 1:2) 0], [0 -1 0])))*acos(dot([alignedBackboneCoordinates(residueIndex*3, 1:2) 0]/norm([alignedBackboneCoordinates(residueIndex*3, 1:2) 0]), [0 -1 0], 2)));
% end