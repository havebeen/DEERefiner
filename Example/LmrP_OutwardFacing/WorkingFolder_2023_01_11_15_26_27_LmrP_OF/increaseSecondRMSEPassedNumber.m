function increaseSecondRMSEPassedNumber(secondRMSEPassedNumberFileName)
            secondRMSEPassedNumber = load(secondRMSEPassedNumberFileName);
            secondRMSEPassedNumber = secondRMSEPassedNumber+1;
            secondRMSEPassedNumberFileIndex = fopen(secondRMSEPassedNumberFileName, 'w');
            fprintf(secondRMSEPassedNumberFileIndex, '%d', secondRMSEPassedNumber);
            fclose(secondRMSEPassedNumberFileIndex);
end