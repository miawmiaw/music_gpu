// Copyright 2017 Chun Shen

#ifndef SRC_FIELD_H_
#define SRC_FIELD_H_

struct Field {
     double *e_rk0;
     double *e_rk1;
     double *e_prev;
     double *e_perm;
     double *rhob_rk0;
     double *rhob_rk1;
     double *rhob_prev;
     double *rhob_perm;
     double **u_rk0;
     double **u_rk1;
     double **u_prev;
     double **u_perm;
     double **dUsup;
     double **Wmunu_rk0;
     double **Wmunu_rk1;
     double **Wmunu_prev;
     double **Wmunu_perm;
     double *pi_b_rk0;
     double *pi_b_rk1;
     double *pi_b_prev;
     double *pi_b_perm;
};

#endif  // SRC_FIELD_H_
  
