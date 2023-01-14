function advancedPrFileCreator(app)
    advanceFileName = char(erase(app.runFileName(1),"runFile_"));
    advanceFileName(end-1:end) = [];
    fid = fopen(strcat(advanceFileName, '_advanced_Pr.dat'), 'w');
    for PrFilesNumber = 1:length(app.DEERRefineRestraintsTable.Data(:, 1))
        fprintf(fid, strcat(advanceFileName, "_", app.DEERRefineRestraintsTable.Data(PrFilesNumber, 2), "_", app.DEERRefineRestraintsTable.Data(PrFilesNumber, 3), "R1.dat\n"));
    end
    fclose(fid);
end