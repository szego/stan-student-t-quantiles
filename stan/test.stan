functions{
  #include student_t_rngs.stan
}

generated quantities {
  real t = student_t_alternate_rng(5, 3, 1);
  real t_lub = student_t_lub_rng(5, 3, 1, 1, 6);
}
