function pdbSaver(new_filnm, formated_structure, states_tag, states_num)
    switch nargin
        case 4
            fid = fopen(new_filnm, 'a');
            fprintf(fid, strcat("MODEL ", num2str(states_num), '\n'));
            fprintf(fid, '%4s  %+5s %-4s%1s%+3s %1s%+4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%+2s%-2s\n', formated_structure');
            fprintf(fid, "ENDMDL\n");
                if states_tag == 2
                    fprintf(fid, "END");
                    fclose(fid);
                end
        case 2
            fid = fopen(new_filnm, 'w');
            fprintf(fid, '%4s  %+5s %-4s%1s%+3s %1s%+4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%+2s%-2s\n', formated_structure');
            fclose(fid);
        otherwise
            fprintf('wrong inputs\n')
    end
end