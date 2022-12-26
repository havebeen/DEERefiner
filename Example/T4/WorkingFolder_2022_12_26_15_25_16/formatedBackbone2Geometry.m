function geometry = formatedBackbone2Geometry(formatedBackbone)
    formatedN  = formatedBackbone((formatedBackbone(:, 3) == " N  "), :);
    formatedCA = formatedBackbone((formatedBackbone(:, 3) == " CA "), :);
    formatedC  = formatedBackbone((formatedBackbone(:, 3) == " C  "), :);
    formatedO  = formatedBackbone((formatedBackbone(:, 3) == " O  "), :);
    
    coorN  = double(formatedN(:, 9:11));
    coorCA = double(formatedCA(:, 9:11));
    coorC  = double(formatedC(:, 9:11));
    coorO  = double(formatedO(:, 9:11));
    
    vecNCA = coorCA-coorN;
    vecCAC = coorC-coorCA;
    vecCO  = coorO-coorC;
    vecCN1 = coorN(2:end, :)-coorC(1:end-1, :);
    
    geometry.length.NCA = sqrt(sum(vecNCA.^2, 2));
    geometry.length.CAC = sqrt(sum(vecCAC.^2, 2));
    geometry.length.CO  = sqrt(sum(vecCO.^2, 2));
    geometry.length.CN1 = sqrt(sum(vecCN1.^2, 2));
    
    normVecCAN = -vecNCA./geometry.length.NCA;
    normVecCAC = vecCAC./geometry.length.CAC;    
    geometry.angle.NCAC = acos(dot(normVecCAN, normVecCAC, 2))/pi*180;
    
    normVecCCA = -vecCAC(1:end-1, :)./geometry.length.CAC(1:end-1);
    normVecCN1 = vecCN1./geometry.length.CN1;
    geometry.angle.CACN1 = acos(dot(normVecCCA, normVecCN1, 2))/pi*180;
    
    normVecN1C = -vecCN1./geometry.length.CN1;
    normVecN1CA1 = vecNCA(2:end, :)./geometry.length.NCA(2:end);
    geometry.angle.CN1CA1 = acos(dot(normVecN1C, normVecN1CA1, 2))/pi*180;
    
    normVecCCA = -vecCAC./geometry.length.CAC;
    normVecCO = vecCO./geometry.length.CO;
    geometry.angle.CACO = acos(dot(normVecCCA, normVecCO, 2))/pi*180;
    
    normVecNCA = vecNCA./geometry.length.NCA;
    normVecCAC = vecCAC./geometry.length.CAC;
    normVecCO = vecCO./geometry.length.CO;
    crossNCAC = cross(-normVecNCA, normVecCAC)./DNorm2(cross(-normVecNCA, normVecCAC), 2);
    crossCACO = cross(-normVecCAC, normVecCO)./DNorm2(cross(-normVecCAC, normVecCO), 2);
    dihedralNCACOValue = acos(dot(crossNCAC, crossCACO, 2))/pi*180;
    dihedralNCACOSignVec = cross(crossNCAC, crossCACO);
    dihedralNCACOSign = sign(dot(dihedralNCACOSignVec, normVecCAC, 2));
    geometry.dihedral.NCACO = dihedralNCACOValue.*dihedralNCACOSign;
    
    normVecNCA = vecNCA(1:end-1, :)./geometry.length.NCA(1:end-1, :);
    normVecCAC = vecCAC(1:end-1, :)./geometry.length.CAC(1:end-1, :);
    normVecCN1 = vecCN1./geometry.length.CN1;
    crossNCAC = cross(-normVecNCA, normVecCAC)./DNorm2(cross(-normVecNCA, normVecCAC), 2);
    crossCACN1 = cross(-normVecCAC, normVecCN1)./DNorm2(cross(-normVecCAC, normVecCN1), 2);
    dihedralNCACN1Value = acos(dot(crossNCAC, crossCACN1, 2))/pi*180;
    dihedralNCACN1SignVec = cross(crossNCAC, crossCACN1);
    dihedralNCACN1Sign = sign(dot(dihedralNCACN1SignVec, normVecCAC, 2));
    geometry.dihedral.NCACN1 = dihedralNCACN1Value.*dihedralNCACN1Sign;
    
    normVecCN1 = vecCN1./geometry.length.CN1;
    normVecN1CA1 = vecNCA(2:end, :)./geometry.length.NCA(2:end);
    normVecCA1C1 = vecCAC(2:end, :)./geometry.length.CAC(2:end);
    crossCN1CA1 = cross(-normVecCN1, normVecN1CA1)./DNorm2(cross(-normVecCN1, normVecN1CA1), 2);
    crossN1CA1C1 = cross(-normVecN1CA1, normVecCA1C1)./DNorm2(cross(-normVecN1CA1, normVecCA1C1), 2);
    dihedralCN1CA1C1Value = acos(dot(crossCN1CA1, crossN1CA1C1, 2))/pi*180;
    dihedralCN1CA1C1SignVec = cross(crossCN1CA1, crossN1CA1C1);
    dihedralCN1CA1C1Sign = sign(dot(dihedralCN1CA1C1SignVec, normVecN1CA1, 2));
    geometry.dihedral.CN1CA1C1 = dihedralCN1CA1C1Value.*dihedralCN1CA1C1Sign;
    
    normVecCAC = vecCAC(1:end-1, :)./geometry.length.CAC(1:end-1, :);
    normVecCN1 = vecCN1./geometry.length.CN1;
    normVecN1CA1 = vecNCA(2:end, :)./geometry.length.NCA(2:end);
    crossCACN1 = cross(-normVecCAC, normVecCN1)./DNorm2(cross(-normVecCAC, normVecCN1), 2);
    crossCN1CA1 = cross(-normVecCN1, normVecN1CA1)./DNorm2(cross(-normVecCN1, normVecN1CA1), 2);
    dihedralCACN1CA1Value = acos(dot(crossCACN1, crossCN1CA1, 2))/pi*180;
    dihedralCACN1CA1SignVec = cross(crossCACN1, crossCN1CA1);
    dihedralCACN1CA1Sign = sign(dot(dihedralCACN1CA1SignVec, normVecCN1, 2));
    geometry.dihedral.CACN1CA1 = dihedralCACN1CA1Value.*dihedralCACN1CA1Sign;
    
end