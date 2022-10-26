function increaseFourthRMSEPassedNumber(fourthRMSEPassedNumberFileName)
            fourthRMSEPassedNumber = load(fourthRMSEPassedNumberFileName);
            fourthRMSEPassedNumber = fourthRMSEPassedNumber+1;
            fourthRMSEPassedNumberFileIndex = fopen(fourthRMSEPassedNumberFileName, 'w');
            fprintf(fourthRMSEPassedNumberFileIndex, '%d', fourthRMSEPassedNumber);
            fclose(fourthRMSEPassedNumberFileIndex);
end