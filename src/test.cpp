
#include <stdio.h>
#include <cmath>

using namespace std;

struct Field {
     double *e_rk0;
     double *e_rk1;
     double *rhob_rk0;
     double *rhob_rk1;
     double **u_rk0;
     double **u_rk1;
};

#pragma acc routine seq
int get_min (int a, int b) {
    if (a > b)
        return (b);
    else
        return(a);
}

#pragma acc routine seq
int get_max (int a, int b) {
    if (a < b)
        return (b);
    else
        return(a);
}

#pragma acc routine seq
void update_grid_array_from_field(
                Field *hydro_fields, int idx, double *grid_array, int rk_flag);
#pragma acc routine seq
void get_qmu_from_grid_array(double tau, double *qi,
                                      double *grid_array);

void initialize_hydro_fields(Field *hydro_fields) {
    int n_cell = 201*201*1;
    hydro_fields->e_rk0 = new double[n_cell];
    hydro_fields->e_rk1 = new double[n_cell];
    hydro_fields->rhob_rk0 = new double[n_cell];
    hydro_fields->rhob_rk1 = new double[n_cell];
    hydro_fields->u_rk0 = new double* [n_cell];
    hydro_fields->u_rk1 = new double* [n_cell];
    for (int i = 0; i < n_cell; i++) {
        hydro_fields->e_rk0[i] = 1.0;
        hydro_fields->e_rk1[i] = 1.0;
        hydro_fields->rhob_rk0[i] = 1.0;
        hydro_fields->rhob_rk1[i] = 1.0;
        hydro_fields->u_rk0[i] = new double[4];
        hydro_fields->u_rk1[i] = new double[4];
    }
}

#pragma acc routine seq
void prepare_qi_array(
        double tau, int rk_flag, Field *hydro_fields, int ieta, int ix, int iy,
        int n_cell_eta, int n_cell_x, int n_cell_y,
        double qi_array[][5], double qi_nbr_x[][5],
        double qi_nbr_y[][5], double qi_nbr_eta[][5],
        double qi_rk0[][5], double grid_array[][5], double *grid_array_temp) {

    double tau_rk;
    if (rk_flag == 0) {
        tau_rk = tau;
    } else {
        tau_rk = tau + 0.1;
    }

    int field_idx;
    int field_ny = 200 + 1;
    int field_nperp = 201*201;
    // first build qi cube n_cell_x*n_cell_x*n_cell_eta
    for (int k = 0; k < n_cell_eta; k++) {
        int idx_ieta = get_min(ieta + k, 0);
        for (int i = 0; i < n_cell_x; i++) {
            int idx_ix = get_min(ix + i, 200);
            for (int j = 0; j < n_cell_y; j++) {
                int idx_iy = get_min(iy + j, 200);
                int idx = j + n_cell_y*i + n_cell_x*n_cell_y*k;
                field_idx = (idx_iy + idx_ix*field_ny + idx_ieta*field_nperp);
                update_grid_array_from_field(hydro_fields, field_idx,
                                             grid_array[idx], rk_flag);
                get_qmu_from_grid_array(tau_rk, qi_array[idx],
                                        grid_array[idx]);
            }
        }
    }

    if (rk_flag == 1) {
        for (int k = 0; k < n_cell_eta; k++) {
            int idx_ieta = get_min(ieta + k, 0);
            for (int i = 0; i < n_cell_x; i++) {
                int idx_ix = get_min(ix + i, 200);
                for (int j = 0; j < n_cell_y; j++) {
                    int idx_iy = get_min(iy + j, 200);
                    int idx = j + n_cell_y*i + n_cell_x*n_cell_y*k;
                    field_idx = (idx_iy + idx_ix*field_ny
                                 + idx_ieta*field_nperp);
                    update_grid_array_from_field(hydro_fields, field_idx,
                                                 grid_array_temp, 0);
                    get_qmu_from_grid_array(tau, qi_rk0[idx],
                                            grid_array_temp);
                }
            }
        }
    }

    // now build neighbouring cells
    // x-direction
    for (int k = 0; k < n_cell_eta; k++) {
        int idx_ieta = get_min(ieta + k, 0);
        for (int i = 0; i < n_cell_y; i++) {
            int idx_iy = get_min(iy + i, 200);
            int idx = 4*i + 4*n_cell_y*k;

            int idx_m_2 = get_max(0, ix - 2);
            int idx_m_1 = get_max(0, ix - 1);
            int idx_p_1 = get_min(ix + n_cell_x, 200);
            int idx_p_2 = get_min(ix + n_cell_x + 1, 200);

            field_idx = (idx_iy + idx_m_2*field_ny + idx_ieta*field_nperp);
            update_grid_array_from_field(hydro_fields, field_idx,
                                         grid_array_temp, rk_flag);
            get_qmu_from_grid_array(tau_rk, qi_nbr_x[idx], grid_array_temp);
            field_idx = (idx_iy + idx_m_1*field_ny + idx_ieta*field_nperp);
            update_grid_array_from_field(hydro_fields, field_idx,
                                         grid_array_temp, rk_flag);
            get_qmu_from_grid_array(tau_rk, qi_nbr_x[idx+1], grid_array_temp);
            field_idx = (idx_iy + idx_p_1*field_ny + idx_ieta*field_nperp);
            update_grid_array_from_field(hydro_fields, field_idx,
                                         grid_array_temp, rk_flag);
            get_qmu_from_grid_array(tau_rk, qi_nbr_x[idx+2], grid_array_temp);
            field_idx = (idx_iy + idx_p_2*field_ny + idx_ieta*field_nperp);
            update_grid_array_from_field(hydro_fields, field_idx,
                                         grid_array_temp, rk_flag);
            get_qmu_from_grid_array(tau_rk, qi_nbr_x[idx+3], grid_array_temp);
        }
    }

    // y-direction
    for (int k = 0; k < n_cell_eta; k++) {
        int idx_ieta = get_min(ieta + k, 1);
        for (int i = 0; i < n_cell_x; i++) {
            int idx_ix = get_min(ix + i, 200);
            int idx = 4*i + 4*n_cell_x*k;

            int idx_m_2 = get_max(0, iy - 2);
            int idx_m_1 = get_max(0, iy - 1);
            int idx_p_1 = get_min(iy + n_cell_y, 200);
            int idx_p_2 = get_min(iy + n_cell_y + 1, 200);

            field_idx = (idx_m_2 + idx_ix*field_ny + idx_ieta*field_nperp);
            update_grid_array_from_field(hydro_fields, field_idx,
                                         grid_array_temp, rk_flag);
            get_qmu_from_grid_array(tau_rk, qi_nbr_y[idx], grid_array_temp);
            field_idx = (idx_m_1 + idx_ix*field_ny + idx_ieta*field_nperp);
            update_grid_array_from_field(hydro_fields, field_idx,
                                         grid_array_temp, rk_flag);
            get_qmu_from_grid_array(tau_rk, qi_nbr_y[idx+1], grid_array_temp);
            field_idx = (idx_p_1 + idx_ix*field_ny + idx_ieta*field_nperp);
            update_grid_array_from_field(hydro_fields, field_idx,
                                         grid_array_temp, rk_flag);
            get_qmu_from_grid_array(tau_rk, qi_nbr_y[idx+2], grid_array_temp);
            field_idx = (idx_p_2 + idx_ix*field_ny + idx_ieta*field_nperp);
            update_grid_array_from_field(hydro_fields, field_idx,
                                         grid_array_temp, rk_flag);
            get_qmu_from_grid_array(tau_rk, qi_nbr_y[idx+3], grid_array_temp);
        }
    }

    // eta-direction
    for (int i = 0; i < n_cell_x; i++) {
        int idx_ix = get_min(ix + i, 200);
        for (int k = 0; k < n_cell_y; k++) {
            int idx_iy = get_min(iy + k, 200);
            int idx = 4*k + 4*n_cell_y*i;

            int idx_m_2 = get_max(0, ieta - 2);
            int idx_m_1 = get_max(0, ieta - 1);
            int idx_p_1 = get_min(ieta + n_cell_eta, 0);
            int idx_p_2 = get_min(ieta + n_cell_eta + 1, 0);

            field_idx = (idx_iy + idx_ix*field_ny + idx_m_2*field_nperp);
            update_grid_array_from_field(hydro_fields, field_idx,
                                         grid_array_temp, rk_flag);
            get_qmu_from_grid_array(tau_rk, qi_nbr_eta[idx], grid_array_temp);
            field_idx = (idx_iy + idx_ix*field_ny + idx_m_1*field_nperp);
            update_grid_array_from_field(hydro_fields, field_idx,
                                         grid_array_temp, rk_flag);
            get_qmu_from_grid_array(tau_rk, qi_nbr_eta[idx+1], grid_array_temp);
            field_idx = (idx_iy + idx_ix*field_ny + idx_p_1*field_nperp);
            update_grid_array_from_field(hydro_fields, field_idx,
                                         grid_array_temp, rk_flag);
            get_qmu_from_grid_array(tau_rk, qi_nbr_eta[idx+2], grid_array_temp);
            field_idx = (idx_iy + idx_ix*field_ny + idx_p_2*field_nperp);
            update_grid_array_from_field(hydro_fields, field_idx,
                                         grid_array_temp, rk_flag);
            get_qmu_from_grid_array(tau_rk, qi_nbr_eta[idx+3], grid_array_temp);
        }
    }
}

void update_grid_array_from_field(
                Field *hydro_fields, int idx, double *grid_array, int rk_flag) {
    if (rk_flag == 0) {
        grid_array[0] = hydro_fields->e_rk0[idx];
        grid_array[4] = hydro_fields->rhob_rk0[idx];
        //for (int i = 1; i < 4; i++) {
        //    grid_array[i] = (hydro_fields->u_rk0[idx][i]
        //                     /hydro_fields->u_rk0[idx][0]);
        //}
    } else {
        grid_array[0] = hydro_fields->e_rk1[idx];
        grid_array[4] = hydro_fields->rhob_rk1[idx];
        //for (int i = 1; i < 4; i++) {
        //    grid_array[i] = (hydro_fields->u_rk1[idx][i]
        //                     /hydro_fields->u_rk1[idx][0]);
        //}
    }
}

void get_qmu_from_grid_array(double tau, double *qi,
                                      double *grid_array) {
    double rhob = grid_array[4];
    double e = grid_array[0];
    //double pressure = get_pressure(e, rhob);
    double pressure = e/3.;
    double gamma = 1./sqrt(1. - grid_array[1]*grid_array[1]
                              - grid_array[2]*grid_array[2]
                              - grid_array[3]*grid_array[3]);
    double gamma_sq = gamma*gamma;
    qi[0] = tau*((e + pressure)*gamma_sq - pressure);
    qi[1] = tau*(e + pressure)*gamma_sq*grid_array[1];
    qi[2] = tau*(e + pressure)*gamma_sq*grid_array[2];
    qi[3] = tau*(e + pressure)*gamma_sq*grid_array[3];
    qi[4] = tau*rhob*gamma;
}


int main() {
    int n_cell_length = 201*201*1;
    Field *hydro_fields = new Field;
    initialize_hydro_fields(hydro_fields);
#pragma acc data copyin (hydro_fields[0:1], \
                         hydro_fields->e_rk0[0:n_cell_length], \
                         hydro_fields->e_rk1[0:n_cell_length], \
                         hydro_fields->rhob_rk0[0:n_cell_length], \
                         hydro_fields->rhob_rk1[0:n_cell_length], \
                         hydro_fields->u_rk0[0:n_cell_length][0:4], \
                         hydro_fields->u_rk1[0:n_cell_length][0:4])
{
    int n_cell_eta = 1;
    int n_cell_x = 1;
    int n_cell_y = 1;
    double tau = 1.0;
    int rk_flag = 0;
    double grid_array[1][5], qi_array[1][5], qi_array_new[1][5], qi_rk0[1][5];
    double qi_nbr_x[4][5], qi_nbr_y[4][5], qi_nbr_eta[4][5];
    double *grid_array_temp = new double[5];
#pragma acc parallel
{
#pragma acc loop private ( qi_array[0:1][0:5],\
                           qi_array_new[0:1][0:5],\
                           qi_rk0[0:1][0:5],\
                           grid_array[0:1][0:5],\
                           qi_nbr_x[0:4][0:5],\
                           qi_nbr_y[0:4][0:5],\
                           qi_nbr_eta[0:4][0:5],\
                           grid_array_temp[0:5])
    for (int ieta = 0; ieta < 1; ieta += n_cell_eta) {
        for (int ix = 0; ix <= 200; ix += n_cell_x) {
            for (int iy = 0; iy <= 200; iy += n_cell_y) {
                prepare_qi_array(tau, rk_flag, hydro_fields, ieta, ix, iy,
                                 n_cell_eta, n_cell_x, n_cell_y, qi_array,
                                 qi_nbr_x, qi_nbr_y, qi_nbr_eta,
                                 qi_rk0, grid_array, grid_array_temp);
            }
        }
    }
}
    delete[] grid_array_temp;
}
    return(0);
}
