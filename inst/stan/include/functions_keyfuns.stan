real int_halfnorm_point(real sigma, real a, real b){
  real s2 = pow(sigma, 2);
  return s2 * ((1 - exp(-b*b / (2*s2))) - (1 - exp(-a*a / (2*s2))));
}

real int_halfnorm_line(real sigma, real a, real b){
  real den = sqrt(2) * sigma;
  return sqrt(pi()/2) * sigma * (erf(b/den) - erf(a/den));
}

real int_halfnorm(real log_sigma, real a, real b, int point){
  real out;
  real sigma = exp(log_sigma);
  if(point){
    out = int_halfnorm_point(sigma, a, b);
  } else{
    out = int_halfnorm_line(sigma, a, b);
  }
  return out;
}

real int_negexp_point(real rate, real a, real b){
  return rate * exp(-a/rate) * (a+rate) -
         rate * exp(-b/rate) * (b+rate);
}

real int_negexp_line(real rate, real a, real b){
  return rate * (exp(-a/rate) - exp(-b/rate));
}

real int_negexp(real log_rate, real a, real b, int point){
  real out;
  real rate = exp(log_rate);
  if(point){
    out = int_negexp_point(rate, a, b);
  } else{
    out = int_negexp_line(rate, a, b);
  }
  return out;
}

real prob_dist(real par1, real par2, int keyfun, real a, real b, int point){
  real out;
  if(keyfun == 0){
    out = int_halfnorm(par1, a, b, point);
  } else if(keyfun == 1){
    out = int_negexp(par1, a, b, point);
  }
  return out;
}
