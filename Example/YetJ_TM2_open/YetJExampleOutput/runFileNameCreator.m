function runFileName = runFileNameCreator(app)
    currentDateAndTime = datestr(now,'YYYY_mm_DD_HH_MM_SS');
    runFileNumber = app.numberOfCPUCores;
    runFileSerial = 1:runFileNumber;
    runFileName = strcat("runFile_", currentDateAndTime, '_', string(runFileSerial));
end