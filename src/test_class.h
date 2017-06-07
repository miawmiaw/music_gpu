#include <stdio.h>

struct Field {
     double *e_rk0;
     double *e_rk1;
     double *rhob_rk0;
     double *rhob_rk1;
     double **u_rk0;
     double **u_rk1;
};

class test{


public:
    test() {};
    ~test() {};

#pragma acc routine seq
    int get_min (int a, int b);

#pragma acc routine seq
    int get_max (int a, int b);

    void initialize_hydro_fields(Field *hydro_fields);

#pragma acc routine seq
    void prepare_qi_array(
        double tau, int rk_flag, Field *hydro_fields, int ieta, int ix, int iy,
        int n_cell_eta, int n_cell_x, int n_cell_y,
        double qi_array[][5], double qi_nbr_x[][5],
        double qi_nbr_y[][5], double qi_nbr_eta[][5],
        double qi_rk0[][5], double grid_array[][5], double *grid_array_temp);

#pragma acc routine seq
    void update_grid_array_from_field(
                Field *hydro_fields, int idx, double *grid_array, int rk_flag);

#pragma acc routine seq
    void get_qmu_from_grid_array(double tau, double *qi,
                                      double *grid_array);

    int run();

};
