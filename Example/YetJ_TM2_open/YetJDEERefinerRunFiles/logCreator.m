function logCreator(app)
    logFileName = char(erase(app.runFileName(1),"runFile_"));
    logFileName(end-1:end) = [];
    
    logFileIndex = fopen(strcat(logFileName, ".log"), 'w');
    
    fprintf(logFileIndex, "advanced_Pr file name\t%s\n", strcat(logFileName, '_advanced_Pr.dat'));
    
    for PrFilesNumber = 1:length(app.DEERRefineRestraintsTable.Data(:, 1))
        residue1 = app.DEERRefineRestraintsTable.Data(PrFilesNumber, 2);
        residue2 = app.DEERRefineRestraintsTable.Data(PrFilesNumber, 3);
        fprintf(logFileIndex, "advanced_Pr file %d\t%s\n", PrFilesNumber, strcat(logFileName, "_", residue1, "_", residue2, "R1.dat"));
        fprintf(logFileIndex, "advanced_Pr original file %d\t%s\n", PrFilesNumber, app.distanceDistributionFullPath(PrFilesNumber));
    end
    
    fprintf(logFileIndex, "Flexbile regions\t%s\n", app.loopRegionString);
    fprintf(logFileIndex, "phi/psi angle\t%f\n", app.phiPsiAngle);
    fprintf(logFileIndex, "MC iteration\t%d\n", app.monteCarloIteration);
    fprintf(logFileIndex, "MC temperature\t%f\n", app.monteCarloTemperature);
    fprintf(logFileIndex, "RMSE\t%f\n", app.RMSE);
    fprintf(logFileIndex, "Maximal clashes\t%d\n", app.maximalClashes);
    fprintf(logFileIndex, "Number of structures\t%d\n", app.numberOfStructure);
    fclose(logFileIndex);
    
end