clear
clc
close all

count = 0;

for i=1:29
    for j=1:10
        for k=1:5
            check = (i - 1)*10*5 + (j - 1)*5 + k;
            count = count + 1;
        end
    end
end

check
count