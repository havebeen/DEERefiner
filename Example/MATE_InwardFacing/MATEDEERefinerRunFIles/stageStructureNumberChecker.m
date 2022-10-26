function stageStructureNumberChecker(app)
    timeStamp = char(erase(app.runFileName(1),"runFile_"));
    timeStamp(end-1:end) = [];
    firstRMSEPassedNumberFileName = strcat(timeStamp, '_first_RMSE_passed_number.dat');
    app.Stage1EditField.Value = load(firstRMSEPassedNumberFileName);
    secondRMSEPassedNumberFileName = strcat(timeStamp, '_second_RMSE_passed_number.dat');
    app.Stage2EditField.Value = load(secondRMSEPassedNumberFileName);
    thirdRMSEPassedNumberFileName = strcat(timeStamp, '_third_RMSE_passed_number.dat');
    app.Stage3EditField.Value = load(thirdRMSEPassedNumberFileName);
    fourthRMSEPassedNumberFileName = strcat(timeStamp, '_fourth_RMSE_passed_number.dat');
    app.Stage4EditField.Value = load(fourthRMSEPassedNumberFileName);
    fifthRMSEPassedNumberFileName = strcat(timeStamp, '_fifth_RMSE_passed_number.dat');
    app.Stage5EditField.Value = load(fifthRMSEPassedNumberFileName);
end