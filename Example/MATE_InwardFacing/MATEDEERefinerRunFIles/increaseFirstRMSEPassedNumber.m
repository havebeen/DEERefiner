function increaseFirstRMSEPassedNumber(firstRMSEPassedNumberFileName)
            firstRMSEPassedNumber = load(firstRMSEPassedNumberFileName);
            firstRMSEPassedNumber = firstRMSEPassedNumber+1;
            firstRMSEPassedNumberFileIndex = fopen(firstRMSEPassedNumberFileName, 'w');
            fprintf(firstRMSEPassedNumberFileIndex, '%d', firstRMSEPassedNumber);
            fclose(firstRMSEPassedNumberFileIndex);
end