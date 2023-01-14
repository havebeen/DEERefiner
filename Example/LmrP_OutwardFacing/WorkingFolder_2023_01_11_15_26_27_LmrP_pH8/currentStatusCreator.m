function currentStatusCreator(app)
    timeStamp = char(erase(app.runFileName(1),"runFile_"));
    timeStamp(end-1:end) = [];
    currentStatusFileName = strcat("currentStatus_", timeStamp, ".dat");
    currentStatusFileIndex = fopen(currentStatusFileName, 'w');
    fprintf(currentStatusFileIndex, "%d", 0);
    fclose(currentStatusFileIndex);
end