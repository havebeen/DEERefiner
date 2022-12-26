function passedStructureNumbers = passedStructureNumbersDetector(app)
    timeStamp = char(erase(app.runFileName(1),"runFile_"));
    timeStamp(end-1:end) = [];firstRMSEPassedNumberFileName = strcat(timeStamp, '_first_RMSE_passed_number.dat');
    passedStructureNumbers(1) = load(firstRMSEPassedNumberFileName);
    secondRMSEPassedNumberFileName = strcat(timeStamp, '_second_RMSE_passed_number.dat');
    passedStructureNumbers(2) = load(secondRMSEPassedNumberFileName);
    thirdRMSEPassedNumberFileName = strcat(timeStamp, '_third_RMSE_passed_number.dat');
    passedStructureNumbers(3) = load(thirdRMSEPassedNumberFileName);
    fourthRMSEPassedNumberFileName = strcat(timeStamp, '_fourth_RMSE_passed_number.dat');
    passedStructureNumbers(4) = load(fourthRMSEPassedNumberFileName);
    fifthRMSEPassedNumberFileName = strcat(timeStamp, '_fifth_RMSE_passed_number.dat');
    passedStructureNumbers(5) = load(fifthRMSEPassedNumberFileName);
end