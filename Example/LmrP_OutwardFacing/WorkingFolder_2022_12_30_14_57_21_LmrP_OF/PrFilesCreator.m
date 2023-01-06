function PrFilesCreator(app)
    advanceFileName = char(erase(app.runFileName(1),"runFile_"));
    advanceFileName(end-1:end) = [];
    for PrFilesNumber = 1:length(app.DEERRefineRestraintsTable.Data(:, 1))
        PrXAxis = reshape(10:90, length(10:90), 1);
        PrYAxis = zeros(length(PrXAxis), 1);
        residue1 = app.DEERRefineRestraintsTable.Data(PrFilesNumber, 2);
        residue2 = app.DEERRefineRestraintsTable.Data(PrFilesNumber, 3);
        targetDistance = round(double(app.DEERRefineRestraintsTable.Data(PrFilesNumber, 4))*10);
        PrYAxis(PrXAxis == targetDistance) = 1;
        PrFileIdex = fopen(strcat(advanceFileName, "_", residue1, "_", residue2, "R1.dat"), 'w');
        fprintf(PrFileIdex, strcat("#", residue1, "\n"));
        fprintf(PrFileIdex, strcat("#", residue2, "\n"));
        fprintf(PrFileIdex, "%d\t%d\n", [PrXAxis PrYAxis]');
        fclose(PrFileIdex);
    end
end