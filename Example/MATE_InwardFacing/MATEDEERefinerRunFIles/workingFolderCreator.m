function workingFolderCreator(app)
    cd(app.initialStructurePath);
    timeStamp = char(erase(app.runFileName(1),"runFile_"));
    timeStamp(end-1:end) = [];
    workingFolderName = strcat("WorkingFolder_", timeStamp);
    mkdir(workingFolderName)
    cd(workingFolderName)
end