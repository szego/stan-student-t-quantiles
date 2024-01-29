vector n_t_q_inner_coeffs(int n, real nu) {
  if (n <= 0)
    reject("n must be positive; found n=", n);
  
  real g = sqrt(nu/2) * exp(lgamma(nu/2) - lgamma((nu+1)/2));
  vector[3] x0 = [
    g,
    ((nu+1)*g^3 - nu*g)/(6*nu),
    ((7*nu^2 + 8*nu + 1)*g^5 - 10*(nu^2 + nu)*g^3 + 3*nu^2*g)/(120*nu^2)
  ]';
  
  if(n < 4)
    return x0[1:n];
  
  vector[n] x = append_row(x0, rep_vector(0, n-3));
  
  for (i in 2:(n-2)) {
    x[i+2] = -(2*i+1)*x[i+1];
    
    for (l in 0:i) {
      for (m in 0:(i-l)) {
        x[i+2] = x[i+2] + ((1+1/nu)*(2*l+1)*(2*m+1) - 2*m*(2*m+1)/nu)*x[i-l-m+1]*x[l+1]*x[m+1];
      }
    }
    for (l in 0:(i-1)) {
      for (m in 0:(i-l-1)) {
        x[i+2] = x[i+2] - (1/nu)*(2*m+1)*x[i-l-m]*x[l+1]*x[m+1];
      }
    }
    
    x[i+2] = x[i+2]/(2*i+3)/(2*i+2);
  }
  
  return x;
}

real n_t_q(real z, real nu) {
  if(abs(z) < 2*sqrt(nu)) {
    vector[11] a = n_t_q_inner_coeffs(11, nu);
    for (i in 0:10)
      a[i+1] *= z^(2*i+1);
    return sum(a);
  } else {
    real w = (1-std_normal_cdf(abs(z)))*nu*sqrt(pi())*exp(lgamma(nu/2) - lgamma((nu+1)/2));
    real sign_z = z > 0 ? 1 : z < 0 ? -1 : 0;
    return sign_z*sqrt(nu)*w^(-1/nu)*(1 - (nu+1)/(2*(nu+2))*w^(2/nu));
  }
}

real student_t_alternate_rng(real nu, real mu, real sigma) {
  return mu + sigma * n_t_q(std_normal_rng(), nu);
}

real student_t_lub_rng(real nu, real mu, real sigma, real lb, real ub) {
  real p_lb = student_t_cdf(lb | nu, mu, sigma);
  real p_ub = student_t_cdf(ub | nu, mu, sigma);
  real u = uniform_rng(p_lb, p_ub);
  real z = inv_Phi(u);
  return mu + sigma * n_t_q(z, nu);
}