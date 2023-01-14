function currentStatusFlag = currentStatusDetector(app)
    timeStamp = char(erase(app.runFileName(1),"runFile_"));
    timeStamp(end-1:end) = [];
    currentStatusFileName = strcat("currentStatus_", timeStamp, ".dat");
    currentStatusFlag = load(currentStatusFileName);
end