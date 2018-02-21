function  f = mask_surrogateK(x, mod1, mod2)
f = zeros(1,2);

f(1) = predictor(x,mod1);
f(2) = predictor(x,mod2);

end
