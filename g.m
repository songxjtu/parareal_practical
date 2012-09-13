function y=g(x)

y=1.6656e-5*exp(25.785*(x-1)./ ...
              x)./(1.6656e-5+ exp(-25.785./x));
end