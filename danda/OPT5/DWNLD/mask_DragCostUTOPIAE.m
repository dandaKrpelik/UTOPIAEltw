function [ f ] = mask_DragCostUTOPIAE(xa)
b = [10,10,100,1000,2000,10000];
x = xa.*b;

f = DragCostUTOPIAE(x);
end