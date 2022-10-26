function flexiblePhiPsiIndex = loopRegionString2FlexiblePhiPsiIndex(app)
    loopRegionStringExecutable = strcat("[", app.loopRegionString, "];");
    loopRegionResidueNumber = eval(loopRegionStringExecutable);
    residueIndex = app.residueIndex;
    [~, flexiblePhiPsiIndex] = intersect(residueIndex,loopRegionResidueNumber);
    
end