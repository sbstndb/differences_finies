function [A, CA, RA] = myldlt_to_delete(SA)
    [n, m] = size(SA);
    A = [];
    CA = [];
    RA = [];
    NNZ = 1;

    RA = [RA NNZ];
    for i = 1 : n
        for j = 1 : m
            if SA(i, j) ~= 0
                A = [A SA(i, j)];
                CA = [CA j];
                NNZ = NNZ + 1;
            end
        end
        RA = [RA NNZ];
    end
end