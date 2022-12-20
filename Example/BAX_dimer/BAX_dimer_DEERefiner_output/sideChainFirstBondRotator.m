function firstBondRotatedSideChain = sideChainFirstBondRotator(formatedSideChain, rotatingAngle)
    residueType = formatedSideChain(1, 5);
    if residueType == "GLY" || residueType == "PRO" || residueType == "ALA"
        firstBondRotatedSideChain = formatedSideChain;
    elseif residueType == "ARG" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        ARGRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD " |...
                     formatedSideChain(:, 3) == " NE " |...
                     formatedSideChain(:, 3) == " CZ " |...
                     formatedSideChain(:, 3) == " NH1" |...
                     formatedSideChain(:, 3) == " NH2", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD " |...
                     formatedSideChain(:, 3) == " NE " |...
                     formatedSideChain(:, 3) == " CZ " |...
                     formatedSideChain(:, 3) == " NH1" |...
                     formatedSideChain(:, 3) == " NH2", :)*ARGRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;
    elseif residueType == "HIS" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        HISRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);  
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD2" |...
                     formatedSideChain(:, 3) == " ND1" |...
                     formatedSideChain(:, 3) == " CE1" |...
                     formatedSideChain(:, 3) == " NE2", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD2" |...
                     formatedSideChain(:, 3) == " ND1" |...
                     formatedSideChain(:, 3) == " CE1" |...
                     formatedSideChain(:, 3) == " NE2", :)*HISRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;
    elseif residueType == "LYS"
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        LYSRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD " |...
                     formatedSideChain(:, 3) == " CE " |...
                     formatedSideChain(:, 3) == " NZ ", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD " |...
                     formatedSideChain(:, 3) == " CE " |...
                     formatedSideChain(:, 3) == " NZ ", :)*LYSRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;
    elseif residueType == "ASP" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        ASPRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " OD1" |...
                     formatedSideChain(:, 3) == " OD2", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " OD1" |...
                     formatedSideChain(:, 3) == " OD2", :)*ASPRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;      
    elseif residueType == "GLU" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        GLURotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD " |...
                     formatedSideChain(:, 3) == " OE1" |...
                     formatedSideChain(:, 3) == " OE2", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD " |...
                     formatedSideChain(:, 3) == " OE1" |...
                     formatedSideChain(:, 3) == " OE2", :)*GLURotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;   
    elseif residueType == "SER"
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        SERRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " OG ", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " OG ", :)*SERRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;   
    elseif residueType == "THR" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        THRRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG2" |...
                     formatedSideChain(:, 3) == " OG1", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG2" |...
                     formatedSideChain(:, 3) == " OG1", :)*THRRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;   
    elseif residueType == "ASN" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        ASNRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " ND2" |...
                     formatedSideChain(:, 3) == " OD1", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " ND2" |...
                     formatedSideChain(:, 3) == " OD1", :)*ASNRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;  
    elseif residueType == "GLN"
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        GLNRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD " |...
                     formatedSideChain(:, 3) == " NE2" |...
                     formatedSideChain(:, 3) == " OE1", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD " |...
                     formatedSideChain(:, 3) == " NE2" |...
                     formatedSideChain(:, 3) == " OE1", :)*GLNRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;  
    elseif residueType == "CYS" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        CYSRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " SG ", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " SG ", :)*CYSRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;  
    elseif residueType == "VAL" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        VALRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CG1" |...
                     formatedSideChain(:, 3) == " CG2", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CG1" |...
                     formatedSideChain(:, 3) == " CG2", :)*VALRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;  
    elseif residueType == "ILE" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        ILERotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG1" |...
                     formatedSideChain(:, 3) == " CG2" |...
                     formatedSideChain(:, 3) == " CD1", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG1" |...
                     formatedSideChain(:, 3) == " CG2" |...
                     formatedSideChain(:, 3) == " CD1", :)*ILERotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;  
    elseif residueType == "LEU"
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        LEURotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD1" |...
                     formatedSideChain(:, 3) == " CD2", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD1" |...
                     formatedSideChain(:, 3) == " CD2", :)*LEURotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;  
    elseif residueType == "MET" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        METRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " SD " |...
                     formatedSideChain(:, 3) == " CE ", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " SD " |...
                     formatedSideChain(:, 3) == " CE ", :)*METRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;  
    elseif residueType == "PHE" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        PHERotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD1" |...
                     formatedSideChain(:, 3) == " CD2" |...
                     formatedSideChain(:, 3) == " CE1" |...
                     formatedSideChain(:, 3) == " CE2" |...
                     formatedSideChain(:, 3) == " CZ ", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD1" |...
                     formatedSideChain(:, 3) == " CD2" |...
                     formatedSideChain(:, 3) == " CE1" |...
                     formatedSideChain(:, 3) == " CE2" |...
                     formatedSideChain(:, 3) == " CZ ", :)*PHERotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;  
    elseif residueType == "TYR"
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        TYRRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD1" |...
                     formatedSideChain(:, 3) == " CD2" |...
                     formatedSideChain(:, 3) == " CE1" |...
                     formatedSideChain(:, 3) == " CE2" |...
                     formatedSideChain(:, 3) == " CZ " |...
                     formatedSideChain(:, 3) == " OH ", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD1" |...
                     formatedSideChain(:, 3) == " CD2" |...
                     formatedSideChain(:, 3) == " CE1" |...
                     formatedSideChain(:, 3) == " CE2" |...
                     formatedSideChain(:, 3) == " CZ " |...
                     formatedSideChain(:, 3) == " OH ", :)*TYRRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;  
    elseif residueType == "TRP" 
        sideChainCoordinates = double(formatedSideChain(:, 9:11));
        sideChainTargetCACoordinates = double(formatedSideChain(formatedSideChain(:, 3) == " CA ", 9:11));
        firstMovedSideChainCoordinates = sideChainCoordinates-sideChainTargetCACoordinates;
        TRPRotatingMatrix = rotatingMatrixGenerator(firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CB ", :), ((rand(1, 1)-0.5)*rotatingAngle)/180*pi);
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD1" |...
                     formatedSideChain(:, 3) == " CD2" |...
                     formatedSideChain(:, 3) == " CE2" |...
                     formatedSideChain(:, 3) == " CE3" |...
                     formatedSideChain(:, 3) == " NE1" |...
                     formatedSideChain(:, 3) == " CZ2" |...
                     formatedSideChain(:, 3) == " CZ3" |...
                     formatedSideChain(:, 3) == " CH2", :) = ...
        firstMovedSideChainCoordinates(formatedSideChain(:, 3) == " CG " |...
                     formatedSideChain(:, 3) == " CD1" |...
                     formatedSideChain(:, 3) == " CD2" |...
                     formatedSideChain(:, 3) == " CE2" |...
                     formatedSideChain(:, 3) == " CE3" |...
                     formatedSideChain(:, 3) == " NE1" |...
                     formatedSideChain(:, 3) == " CZ2" |...
                     formatedSideChain(:, 3) == " CZ3" |...
                     formatedSideChain(:, 3) == " CH2", :)*TRPRotatingMatrix;
        secondMovedSideChainCoordinates = firstMovedSideChainCoordinates+sideChainTargetCACoordinates;
        firstBondRotatedSideChainTmp = formatedSideChain;
        firstBondRotatedSideChainTmp(:, 9:11, :) = string(secondMovedSideChainCoordinates);
        firstBondRotatedSideChain = firstBondRotatedSideChainTmp;  
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    end
end