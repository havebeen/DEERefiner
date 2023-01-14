function RMSEPassedNumberFilesCreator(app)
    timeStamp = char(erase(app.runFileName(1),"runFile_"));
    timeStamp(end-1:end) = [];
    
    firstRMSEPassedNumberFileName = strcat(timeStamp, '_first_RMSE_passed_number.dat');
    firstRMSEPassedNumberFileIndex = fopen(firstRMSEPassedNumberFileName, 'w');
    fprintf(firstRMSEPassedNumberFileIndex, '%d', 0);
    fclose(firstRMSEPassedNumberFileIndex);
    
    secondRMSEPassedNumberFileName = strcat(timeStamp, '_second_RMSE_passed_number.dat');
    secondRMSEPassedNumberFileIndex = fopen(secondRMSEPassedNumberFileName, 'w');
    fprintf(secondRMSEPassedNumberFileIndex, '%d', 0);
    fclose(secondRMSEPassedNumberFileIndex);
    
    thirdRMSEPassedNumberFileName = strcat(timeStamp, '_third_RMSE_passed_number.dat');
    thirdRMSEPassedNumberFileIndex = fopen(thirdRMSEPassedNumberFileName, 'w');
    fprintf(thirdRMSEPassedNumberFileIndex, '%d', 0);
    fclose(thirdRMSEPassedNumberFileIndex);
    
    fourthRMSEPassedNumberFileName = strcat(timeStamp, '_fourth_RMSE_passed_number.dat');
    fourthRMSEPassedNumberFileIndex = fopen(fourthRMSEPassedNumberFileName, 'w');
    fprintf(fourthRMSEPassedNumberFileIndex, '%d', 0);
    fclose(fourthRMSEPassedNumberFileIndex);
    
    fifthRMSEPassedNumberFileName = strcat(timeStamp, '_fifth_RMSE_passed_number.dat');
    fifthRMSEPassedNumberFileIndex = fopen(fifthRMSEPassedNumberFileName, 'w');
    fprintf(fifthRMSEPassedNumberFileIndex, '%d', 0);
    fclose(fifthRMSEPassedNumberFileIndex);
    
end