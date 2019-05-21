clc
close all

dim1 = Output(1, 1);
dim2 = Output(2, 1);
dim3 = Output(3, 1);
dim4 = Output(4, 1);

dV = zeros(dim1, dim2, dim3, dim4);

count = 1;
for i=1:dim1
    for j=1:dim2
        for k=1:dim3
            for l=1:dim4
                dV(i, j, k, l) = Output(5, count);
                count = count+1;
            end
        end
    end
end

phi = 10:10:200;
dV1 = dV(:, 1, 1, 1)

figure
plot(phi, dV1)