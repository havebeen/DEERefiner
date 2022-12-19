function DEERefineMFileChecker(app)
    DEERFileMatrix = 'DEERefineFilesName.mat';
    load(DEERFileMatrix);
    for checkingFIleNumber = 1:length(DEERefineFilesName)
        DEERefineFileExistence(checkingFIleNumber) = exist(DEERefineFilesName(checkingFIleNumber), 'file');
    end
    DEERefineFileAllExistence = all(DEERefineFileExistence);
    if DEERefineFileAllExistence == 1
        app.WaitingformfilesdelectedCheckBox.Value = 1;
        app.WaitingformfilesdelectedCheckBox.Text = "All .m files delected";
        app.WaitingformfilesdelectedCheckBox.FontColor = [0.47 0.67 0.19];
    else
        return
    end
end