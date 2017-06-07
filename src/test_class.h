#include <stdio.h>

struct Field {
     double *e_rk0;
     double *e_rk1;
     double *e_prev;
     double *rhob_rk0;
     double *rhob_rk1;
     double *rhob_prev;
     double **u_rk0;
     double **u_rk1;
     double **u_prev;
     double **dUsup;
     double **Wmunu_rk0;
     double **Wmunu_rk1;
     double **Wmunu_prev;
     double *pi_b_rk0;
     double *pi_b_rk1;
     double *pi_b_prev;
};

class test{

public:
    test() {};
    ~test() {};

#pragma acc routine seq
    int get_min (int a, int b);

#pragma acc routine seq
    int get_max (int a, int b);

#pragma acc routine seq
    double get_pressure(double e_local, double rhob);

#pragma acc routine seq
    void revert_grid(double *grid_array, double *grid_prev);

#pragma acc routine seq
    int ReconstIt_velocity_Newton(
        double *grid_array, double tau, double *uq, double *grid_array_p);

#pragma acc routine seq
    double MaxSpeed(double tau, int direc, double *grid_array);

    void initialize_hydro_fields(Field *hydro_fields);

#pragma acc routine seq
    void prepare_qi_array(
        double tau, int rk_flag, Field *hydro_fields, int ieta, int ix, int iy,
        int n_cell_eta, int n_cell_x, int n_cell_y,
        double qi_array[][5], double qi_nbr_x[][5],
        double qi_nbr_y[][5], double qi_nbr_eta[][5],
        double qi_rk0[][5], double grid_array[][5], double *grid_array_temp);

#pragma acc routine seq
    void MakeDeltaQI(double tau, double qi_array[][5], double qi_nbr_x[][5],
                     double qi_nbr_y[][5], double qi_nbr_eta[][5],
                     int n_cell_eta, int n_cell_x, int n_cell_y,
                     double qi_array_new[][5], double grid_array[][5],
                     double *rhs, double *qiphL, double *qiphR,
                     double *qimhL, double *qimhR,
                     double *grid_array_hL, double *grid_array_hR);

#pragma acc routine seq
    void update_grid_array_from_field(
                Field *hydro_fields, int idx, double *grid_array, int rk_flag);

#pragma acc routine seq
    void get_qmu_from_grid_array(double tau, double *qi,
                                      double *grid_array);

#pragma acc routine seq
    double minmod_dx(double up1, double u, double um1);

#pragma acc routine seq
    double get_TJb_new(double *grid_array, int mu, int nu);

#pragma acc routine seq
    double reconst_u0_df(double u0, double T00, double K00, double M,
                              double J0);
#pragma acc routine seq
    double reconst_u0_f_Newton(double u0, double T00, double K00,
                                    double M, double J0);
#pragma acc routine seq
    double reconst_u0_f(double u0, double T00, double K00, double M,
                             double J0);

#pragma acc routine seq
    double reconst_velocity_df(double v, double T00, double M,
                                    double J0);

#pragma acc routine seq
    double reconst_velocity_f_Newton(double v, double T00, double M,
                                          double J0);

#pragma acc routine seq
    double reconst_velocity_f(double v, double T00, double M,
                                   double J0);
#pragma acc routine seq
    int FirstRKStepT(double tau, int rk_flag,
                          double qi_array[][5], double qi_nbr_x[][5],
                          double qi_nbr_y[][5], double qi_nbr_eta[][5],
                          int n_cell_eta, int n_cell_x, int n_cell_y,
                          double vis_array[][19], double vis_nbr_tau[][19],
                          double vis_nbr_x[][19], double vis_nbr_y[][19],
                          double vis_nbr_eta[][19], double qi_rk0[][5],
                          double qi_array_new[][5], double grid_array[][5],
                          double *rhs, double *qiphL, double *qiphR,
                          double *qimhL, double *qimhR,
                          double *grid_array_hL, double *grid_array_hR);

#pragma acc routine seq
    double p_e_func(double e_local, double rhob);

#pragma acc routine seq
    double p_rho_func(double e_local, double rhob);

#pragma acc routine seq
    double get_mu(double e_local, double rhob);

#pragma acc routine seq
    double get_cs2(double e_local, double rhob);

#pragma acc routine seq
    void update_grid_cell(double grid_array[][5], Field *hydro_fields, int rk_flag, int ieta, int ix, int iy, int n_cell_eta, int n_cell_x, int n_cell_y);

#pragma acc routine seq
    void update_grid_array_to_hydro_fields(
            double *grid_array, Field *hydro_fields, int idx, int rk_flag);

    int run();

};
