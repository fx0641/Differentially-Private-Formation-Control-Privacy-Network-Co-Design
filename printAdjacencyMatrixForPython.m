function [string1,string2] = printAdjacencyMatrixForPython(L,epsilon)
    adj = abs(L-diag(diag(L)));
    [m,n] = size(L);
    string1 = '[';
    string2 = '[';
    for i = 1:m
        rowstring = '[';
        for j = 1:n
            if j < n
                val = num2str(adj(i,j));
                valncomma = append(val,',');
                rowstring = append(rowstring,valncomma);
            else
                val = num2str(adj(i,j));
                valncomma = append(val,']');
                rowstring = append(rowstring,valncomma);
            end
        end
        
        if i < m
            string1 = append(string1,rowstring);
            string1 = append(string1,',');
            string2 = append(string2,num2str(epsilon(i)));
            string2 = append(string2,',');
        else
            string1 = append(string1,rowstring);
            string1 = append(string1,']');
            string2 = append(string2,num2str(epsilon(i)));
            string2 = append(string2,']');
    end
end

