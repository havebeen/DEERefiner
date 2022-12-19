function increaseFifthRMSEPassedNumber(fifthRMSEPassedNumberFileName)
            fifthRMSEPassedNumber = load(fifthRMSEPassedNumberFileName);
            fifthRMSEPassedNumber = fifthRMSEPassedNumber+1;
            fifthRMSEPassedNumberFileIndex = fopen(fifthRMSEPassedNumberFileName, 'w');
            fprintf(fifthRMSEPassedNumberFileIndex, '%d', fifthRMSEPassedNumber);
            fclose(fifthRMSEPassedNumberFileIndex);
end