clear
clc
close all

out = csvread("Output.csv");

dim1 = out(1)
dim2 = out(2)
dim3 = out(3)
dim4 = out(4)

nums = out(5:end);
dV = zeros(dim1, dim2, dim3, dim4);

count = 1;
for i=1:dim1
    for j=1:dim2
        for k=1:dim3
            for l=1:dim4
                ans = nums(count);
                
                if ans == 0
                    a = 8;
                end
                
                dV(i, j, k, l) = ans;
                count = count+1;
            end
        end
    end
end

mindV = min(min(min(min(dV))));

for i=1:dim1
    for j=1:dim2  
        for k=1:dim3
            for l=1:dim4
                test = dV(i, j, k, l);
                if (test == mindV)
                    cords = [i, j, k, l];
                    break
                end
            end
        end
    end
end

