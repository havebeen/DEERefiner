function runFilesCreator(app)
    fileNumbers = length(app.runFileName);
    templateFileName = "appTemplate.m";
    for fileNameNumber = 1:fileNumbers
        fileName(fileNameNumber) = strcat(app.runFileName(fileNameNumber), ".m");
        copyfile(templateFileName, fileName(fileNameNumber))
    end
    
end