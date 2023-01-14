function increaseThirdRMSEPassedNumber(thirdRMSEPassedNumberFileName)
            thirdRMSEPassedNumber = load(thirdRMSEPassedNumberFileName);
            thirdRMSEPassedNumber = thirdRMSEPassedNumber+1;
            thirdRMSEPassedNumberFileIndex = fopen(thirdRMSEPassedNumberFileName, 'w');
            fprintf(thirdRMSEPassedNumberFileIndex, '%d', thirdRMSEPassedNumber);
            fclose(thirdRMSEPassedNumberFileIndex);
end