function m = quadratic_model(s,f0,df,hf)
m =f0+s*df+.5*s*hf*s';

