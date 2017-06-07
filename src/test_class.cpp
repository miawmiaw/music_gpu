#include "test_class.h"
#include<iostream>
#include<cmath>

using namespace std;

#define maxi(a, b) ((a) > (b) ? (a) : (b))
#define mini(a, b) ((a) < (b) ? (a) : (b))

int test::get_min (int a, int b) {
    if (a > b)
        return (b);
    else
        return(a);
}

int test::get_max (int a, int b) {
    if (a < b)
        return (b);
    else
        return(a);
};

void test::initialize_hydro_fields(Field *hydro_fields) {
    int n_cell = 201*201*1;
    hydro_fields->e_rk0 = new double[n_cell];
    hydro_fields->e_rk1 = new double[n_cell];
    hydro_fields->e_prev = new double[n_cell];
    hydro_fields->rhob_rk0 = new double[n_cell];
    hydro_fields->rhob_rk1 = new double[n_cell];
    hydro_fields->rhob_prev = new double[n_cell];
    hydro_fields->u_rk0 = new double* [n_cell];
    hydro_fields->u_rk1 = new double* [n_cell];
    hydro_fields->u_prev = new double* [n_cell];
    hydro_fields->dUsup = new double* [n_cell];
    hydro_fields->Wmunu_rk0 = new double* [n_cell];
    hydro_fields->Wmunu_rk1 = new double* [n_cell];
    hydro_fields->Wmunu_prev = new double* [n_cell];
    hydro_fields->pi_b_rk0 = new double [n_cell];
    hydro_fields->pi_b_rk1 = new double [n_cell];
    hydro_fields->pi_b_prev = new double [n_cell];
    for (int i = 0; i < n_cell; i++) {
        hydro_fields->e_rk0[i] = 1.0;
        hydro_fields->e_rk1[i] = 1.0;
        hydro_fields->e_prev[i] = 1.0;
        hydro_fields->rhob_rk0[i] = 1.0;
        hydro_fields->rhob_rk1[i] = 1.0;
        hydro_fields->rhob_prev[i] = 1.0;
        hydro_fields->u_rk0[i] = new double[4];
        hydro_fields->u_rk1[i] = new double[4];
        hydro_fields->u_prev[i] = new double[4];
        hydro_fields->dUsup[i] = new double [20];
        hydro_fields->Wmunu_rk0[i] = new double[4];
        hydro_fields->Wmunu_rk1[i] = new double[4];
        hydro_fields->Wmunu_prev[i] = new double[4];
    }
}

int test::run() {
    int n_cell_length = 201*201*1;
    Field *hydro_fields = new Field;
    initialize_hydro_fields(hydro_fields);
#pragma acc data copyin (hydro_fields[0:1], \
                         hydro_fields->e_rk0[0:n_cell_length], \
                         hydro_fields->e_rk1[0:n_cell_length], \
                         hydro_fields->e_prev[0:n_cell_length], \
                         hydro_fields->rhob_rk0[0:n_cell_length], \
                         hydro_fields->rhob_rk1[0:n_cell_length], \
                         hydro_fields->rhob_prev[0:n_cell_length], \
                         hydro_fields->u_rk0[0:n_cell_length][0:4], \
                         hydro_fields->u_rk1[0:n_cell_length][0:4], \
                         hydro_fields->u_prev[0:n_cell_length][0:4], \
                         hydro_fields->dUsup[0:n_cell_length][0:20], \
                         hydro_fields->Wmunu_rk0[0:n_cell_length][0:14], \
                         hydro_fields->Wmunu_rk1[0:n_cell_length][0:14], \
                         hydro_fields->Wmunu_prev[0:n_cell_length][0:14], \
                         hydro_fields->pi_b_rk0[0:n_cell_length], \
                         hydro_fields->pi_b_rk1[0:n_cell_length], \
                         hydro_fields->pi_b_prev[0:n_cell_length], \
                         )
{
    double tau0 = 1.;
    double dt = 0.01;
    int itmax = 200;
    for (int it = 0; it <= itmax; it++) {
        double tau = tau0 + dt*it;

        int n_cell_eta = 1;
        int n_cell_x = 1;
        int n_cell_y = 1;
        int rk_flag = 0;
        double grid_array[1][5], qi_array[1][5], qi_array_new[1][5], qi_rk0[1][5];
        double qi_nbr_x[4][5], qi_nbr_y[4][5], qi_nbr_eta[4][5];
        double *grid_array_temp = new double[5];
        double *rhs = new double[5];
        double *qiphL = new double[5];
        double *qiphR = new double[5];
        double *qimhL = new double[5];
        double *qimhR = new double[5];
        double *grid_array_hL = new double[5];
        double *grid_array_hR = new double[5];
        double vis_array[1][19], vis_array_new[1][19], vis_nbr_tau[1][19];
        double velocity_array[1][20];
        double vis_nbr_x[4][19], vis_nbr_y[4][19], vis_nbr_eta[4][19];
#pragma acc parallel
{
#pragma acc loop private ( qi_array[0:1][0:5],\
                           qi_array_new[0:1][0:5],\
                           qi_rk0[0:1][0:5],\
                           grid_array[0:1][0:5],\
                           qi_nbr_x[0:4][0:5],\
                           qi_nbr_y[0:4][0:5],\
                           qi_nbr_eta[0:4][0:5],\
                           vis_array[0:1][0:19],\
                           vis_array_new[0:1][0:19],\
                           vis_nbr_tau[0:1][0:19],\
                           velocity_array[0:1][0:20],\
                           vis_nbr_x[0:4][0:19],\
                           vis_nbr_y[0:4][0:19],\
                           vis_nbr_eta[0:4][0:19], \
                           grid_array_temp[0:5], \
                           rhs[0:5], qiphL[0:5], qiphR[0:5], \
                           qimhL[0:5], qimhR[0:5],\
                           grid_array_hL[0:5], grid_array_hR[0:5])
    for (int ieta = 0; ieta < 1; ieta += n_cell_eta) {
        for (int ix = 0; ix <= 200; ix += n_cell_x) {
            for (int iy = 0; iy <= 200; iy += n_cell_y) {
                prepare_qi_array(tau, rk_flag, hydro_fields, ieta, ix, iy,
                                 n_cell_eta, n_cell_x, n_cell_y, qi_array,
                                 qi_nbr_x, qi_nbr_y, qi_nbr_eta,
                                 qi_rk0, grid_array, grid_array_temp);
                FirstRKStepT(tau, rk_flag,
                             qi_array, qi_nbr_x, qi_nbr_y, qi_nbr_eta,
                             n_cell_eta, n_cell_x, n_cell_y,
                             vis_array, vis_nbr_tau,
                             vis_nbr_x, vis_nbr_y, vis_nbr_eta,
                             qi_rk0, qi_array_new, grid_array,
                             rhs, qiphL, qiphR, qimhL, qimhR,
                             grid_array_hL, grid_array_hR);
            }
        }
    }
}
#pragma acc parallel
{
#pragma acc loop private ( qi_array[0:1][0:5],\
                           qi_array_new[0:1][0:5],\
                           qi_rk0[0:1][0:5],\
                           grid_array[0:1][0:5],\
                           qi_nbr_x[0:4][0:5],\
                           qi_nbr_y[0:4][0:5],\
                           qi_nbr_eta[0:4][0:5],\
                           vis_array[0:1][0:19],\
                           vis_array_new[0:1][0:19],\
                           vis_nbr_tau[0:1][0:19],\
                           velocity_array[0:1][0:20],\
                           vis_nbr_x[0:4][0:19],\
                           vis_nbr_y[0:4][0:19],\
                           vis_nbr_eta[0:4][0:19], \
                           grid_array_temp[0:5], \
                           rhs[0:5], qiphL[0:5], qiphR[0:5], \
                           qimhL[0:5], qimhR[0:5],\
                           grid_array_hL[0:5], grid_array_hR[0:5])
    for (int ieta = 0; ieta < 1; ieta += n_cell_eta) {
        for (int ix = 0; ix <= 200; ix += n_cell_x) {
            for (int iy = 0; iy <= 200; iy += n_cell_y) {
                update_grid_cell(grid_array, hydro_fields, rk_flag,
                                 ieta, ix, iy, n_cell_eta, n_cell_x, n_cell_y);
            }
        }
    }
}
    delete[] grid_array_temp;
    delete[] rhs;
    delete[] qiphL;
    delete[] qiphR;
    delete[] qimhL;
    delete[] qimhR;
    delete[] grid_array_hL;
    delete[] grid_array_hR;
    
    cout << "Done time step " << it << "/" << itmax << ". tau = "
         << tau << " fm/c" << endl;
    }
}
    return(0);
};

/* %%%%%%%%%%%%%%%%%%%%%% First steps begins here %%%%%%%%%%%%%%%%%% */
int test::FirstRKStepT(double tau, int rk_flag,
                          double qi_array[][5], double qi_nbr_x[][5],
                          double qi_nbr_y[][5], double qi_nbr_eta[][5],
                          int n_cell_eta, int n_cell_x, int n_cell_y,
                          double vis_array[][19], double vis_nbr_tau[][19],
                          double vis_nbr_x[][19], double vis_nbr_y[][19],
                          double vis_nbr_eta[][19], double qi_rk0[][5],
                          double qi_array_new[][5], double grid_array[][5],
                          double *rhs, double *qiphL, double *qiphR,
                          double *qimhL, double *qimhR,
                          double *grid_array_hL, double *grid_array_hR) {

    // this advances the ideal part
    double tau_next = tau + 0.1;
    double tau_rk;
    if (rk_flag == 0) {
        tau_rk = tau;
    } else {
        tau_rk = tau_next;
    }

    // Solve partial_a T^{a mu} = -partial_a W^{a mu}
    // Update T^{mu nu}

    // MakeDelatQI gets
    //   qi = q0 if rk_flag = 0 or
    //   qi = q0 + k1 if rk_flag = 1
    // rhs[alpha] is what MakeDeltaQI outputs. 
    // It is the spatial derivative part of partial_a T^{a mu}
    // (including geometric terms)
    MakeDeltaQI(tau_rk, qi_array, qi_nbr_x, qi_nbr_y, qi_nbr_eta,
                n_cell_eta, n_cell_x, n_cell_y, qi_array_new, grid_array,
                rhs, qiphL, qiphR, qimhL, qimhR, grid_array_hL, grid_array_hR);

    // now MakeWSource returns partial_a W^{a mu}
    // (including geometric terms) 
    //MakeWSource(tau_rk, qi_array, n_cell_eta, n_cell_x, n_cell_y,
    //            vis_array, vis_nbr_tau, vis_nbr_x, vis_nbr_y,
    //            vis_nbr_eta, qi_array_new, DATA);
    
    if (rk_flag == 1) {
        // if rk_flag == 1, we now have q0 + k1 + k2. 
        // So add q0 and multiply by 1/2
        for (int k = 0; k < n_cell_eta; k++) {
            for (int i = 0; i < n_cell_x; i++) {
                for (int j = 0; j < n_cell_y; j++) {
                    int idx = j + i*n_cell_y + k*n_cell_y*n_cell_x;
                    for (int alpha = 0; alpha < 5; alpha++) {
                        qi_array_new[idx][alpha] += qi_rk0[idx][alpha];
                        qi_array_new[idx][alpha] *= 0.5;
                    }
                }
            }
        }
    }

    for (int k = 0; k < n_cell_eta; k++) {
        for (int i = 0; i < n_cell_x; i++) {
            for (int j = 0; j < n_cell_y; j++) {
                int idx = j + i*n_cell_y + k*n_cell_x*n_cell_y;
                ReconstIt_velocity_Newton(grid_array[idx], tau_next,
                                          qi_array_new[idx],
                                          grid_array[idx]);
            }
        }
    }
    return(0);
}

void test::get_qmu_from_grid_array(double tau, double *qi,
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
};

void test::update_grid_array_from_field(
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
};

void test::prepare_qi_array(
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

void test::MakeDeltaQI(double tau, double qi_array[][5], double qi_nbr_x[][5],
                          double qi_nbr_y[][5], double qi_nbr_eta[][5],
                          int n_cell_eta, int n_cell_x, int n_cell_y,
                          double qi_array_new[][5], double grid_array[][5],
                          double *rhs, double *qiphL, double *qiphR,
                          double *qimhL, double *qimhR,
                          double *grid_array_hL, double *grid_array_hR) {
    /* \partial_tau (tau Ttautau) + \partial_eta Tetatau 
            + \partial_x (tau Txtau) + \partial_y (tau Tytau) + Tetaeta = 0 */
    /* \partial_tau (tau Ttaueta) + \partial_eta Teteta 
            + \partial_x (tau Txeta) + \partial_y (tau Txeta) + Tetatau = 0 */
    /* \partial_tau (tau Txtau) + \partial_eta Tetax + \partial_x tau T_xx
            + \partial_y tau Tyx = 0 */
    
    //double check0 = 0.0;
    //for (int i = 0; i < 5; i++) {
    //    check0 += fabs(qi_array[0][i] - qi[i]);
    //}
    //cout << check0 << endl;

    // tau*Tmu0
    //double rhs[5];
    for (int alpha = 0; alpha < 5; alpha++) {
        rhs[alpha] = 0.0;
    }

    //double *qiphL = new double[5];
    //double *qiphR = new double[5];
    //double *qimhL = new double[5];
    //double *qimhR = new double[5];
    //
    //double *grid_array_hL = new double[5];
    //double *grid_array_hR = new double[5];
    
    for (int k = 0; k < n_cell_eta; k++) {
        for (int i = 0; i < n_cell_x; i++) {
            for (int j = 0; j < n_cell_y; j++) {
                int idx = j + i*n_cell_y + k*n_cell_x*n_cell_y;

                // implement Kurganov-Tadmor scheme
                // here computes the half way T^\tau\mu currents
                // x-direction
                int direc = 1;
                double tau_fac = tau;
                for (int alpha = 0; alpha < 5; alpha++) {
                    double gp = qi_array[idx][alpha];
                    double gphL = qi_array[idx][alpha];
                    double gmhR = qi_array[idx][alpha];
                    
                    double gphR, gmhL, gphR2, gmhL2;
                    if (i + 1 < n_cell_x) {
                        int idx_p_1 = j + (i+1)*n_cell_y + k*n_cell_x*n_cell_y;
                        gphR = qi_array[idx_p_1][alpha];
                    } else {
                        int idx_p_1 = 4*j + k*4*n_cell_y + 2;
                        gphR = qi_nbr_x[idx_p_1][alpha];
                    }
                    if (i - 1 >= 0) {
                        int idx_m_1 = j + (i-1)*n_cell_y + k*n_cell_x*n_cell_y;
                        gmhL = qi_array[idx_m_1][alpha];
                    } else {
                        int idx_m_1 = 4*j + k*4*n_cell_y + 1;
                        gmhL = qi_nbr_x[idx_m_1][alpha];
                    }
                    if (i + 2 < n_cell_x) {
                        int idx_p_2 = j + (i+2)*n_cell_y + k*n_cell_x*n_cell_y;
                        gphR2 = qi_array[idx_p_2][alpha];
                    } else {
                        int idx_p_2 = 4*j + k*4*n_cell_y + 4 + i - n_cell_x;
                        gphR2 = qi_nbr_x[idx_p_2][alpha];
                    }
                    if (i - 2 >= 0) {
                        int idx_m_2 = j + (i-2)*n_cell_y + k*n_cell_x*n_cell_y;
                        gmhL2 = qi_array[idx_m_2][alpha];
                    } else {
                        int idx_m_2 = 4*j + k*4*n_cell_y + i;
                        gmhL2 = qi_nbr_x[idx_m_2][alpha];
                    }

                    double fphL = 0.5*minmod_dx(gphR, gp, gmhL);
                    double fphR = -0.5*minmod_dx(gphR2, gphR, gp);
                    double fmhL = 0.5*minmod_dx(gp, gmhL, gmhL2);
                    double fmhR = -0.5*minmod_dx(gphR, gp, gmhL);
                    qiphL[alpha] = gphL + fphL;
                    qiphR[alpha] = gphR + fphR;
                    qimhL[alpha] = gmhL + fmhL;
                    qimhR[alpha] = gmhR + fmhR;
                }
                // for each direction, reconstruct half-way cells
                // reconstruct e, rhob, and u[4] for half way cells
                int flag = ReconstIt_velocity_Newton(
                                grid_array_hL, tau, qiphL, grid_array[idx]);
                double aiphL = MaxSpeed(tau, direc, grid_array_hL);

                flag *= ReconstIt_velocity_Newton(
                                grid_array_hR, tau, qiphR, grid_array[idx]);
                double aiphR = MaxSpeed(tau, direc, grid_array_hR);
                double aiph = maxi(aiphL, aiphR);
                for (int alpha = 0; alpha < 5; alpha++) {
                    double FiphL = tau_fac*get_TJb_new(grid_array_hL,
                                                       alpha, direc);
                    double FiphR = tau_fac*get_TJb_new(grid_array_hR,
                                                        alpha, direc);
                    // KT: H_{j+1/2} = (f(u^+_{j+1/2}) + f(u^-_{j+1/2})/2
                    //              - a_{j+1/2}(u_{j+1/2}^+ - u^-_{j+1/2})/2
                    double Fiph = 0.5*((FiphL + FiphR)
                                        - aiph*(qiphR[alpha] - qiphL[alpha]));

                    rhs[alpha] = -Fiph/0.1*0.01;
                }

                flag *= ReconstIt_velocity_Newton(grid_array_hL, tau,
                                                     qimhL, grid_array[idx]);
                double aimhL = MaxSpeed(tau, direc, grid_array_hL);

                flag *= ReconstIt_velocity_Newton(grid_array_hR, tau,
                                                     qimhR, grid_array[idx]);
                double aimhR = MaxSpeed(tau, direc, grid_array_hR);
                double aimh = maxi(aimhL, aimhR);

                for (int alpha = 0; alpha < 5; alpha++) {
                    double FimhL = tau_fac*get_TJb_new(grid_array_hL,
                                                       alpha, direc);
                    double FimhR = tau_fac*get_TJb_new(grid_array_hR,
                                                       alpha, direc);
                    // KT: H_{j+1/2} = (f(u^+_{j+1/2}) + f(u^-_{j+1/2})/2
                    //              - a_{j+1/2}(u_{j+1/2}^+ - u^-_{j+1/2})/2
                    double Fimh = 0.5*((FimhL + FimhR)
                                        - aimh*(qimhR[alpha] - qimhL[alpha]));
                    rhs[alpha] += Fimh/0.1*0.01;
                }
                //cout << "x-direction" << endl;
                
                // y-direction
                direc = 2;
                tau_fac = tau;
                for (int alpha = 0; alpha < 5; alpha++) {
                    double gp = qi_array[idx][alpha];
                    double gphL = qi_array[idx][alpha];
                    double gmhR = qi_array[idx][alpha];

                    double gphR, gmhL, gphR2, gmhL2;
                    if (j + 1 < n_cell_y) {
                        int idx_p_1 = j + 1 + i*n_cell_y + k*n_cell_x*n_cell_y;
                        gphR = qi_array[idx_p_1][alpha];
                    } else {
                        int idx_p_1 = 4*i + 4*k*n_cell_x + 2;
                        gphR = qi_nbr_y[idx_p_1][alpha];
                    }
                    if (j - 1 >= 0) {
                        int idx_m_1 = j - 1 + i*n_cell_y + k*n_cell_x*n_cell_y;
                        gmhL = qi_array[idx_m_1][alpha];
                    } else {
                        int idx_m_1 = 4*i + 4*k*n_cell_x + 1;
                        gmhL = qi_nbr_y[idx_m_1][alpha];
                    }
                    if (j + 2 < n_cell_y) {
                        int idx_p_2 = j + 2 + i*n_cell_y + k*n_cell_x*n_cell_y;
                        gphR2 = qi_array[idx_p_2][alpha];
                    } else {
                        int idx_p_2 = 4*i + 4*k*n_cell_x + 4 + j - n_cell_y;
                        gphR2 = qi_nbr_y[idx_p_2][alpha];
                    }
                    if (j - 2 >= 0) {
                        int idx_m_2 = j - 2 + i*n_cell_y + k*n_cell_x*n_cell_y;
                        gmhL2 = qi_array[idx_m_2][alpha];
                    } else {
                        int idx_m_2 = 4*i + 4*k*n_cell_x + j;
                        gmhL2 = qi_nbr_y[idx_m_2][alpha];
                    }

                    double fphL = 0.5*minmod_dx(gphR, gp, gmhL);
                    double fphR = -0.5*minmod_dx(gphR2, gphR, gp);
                    double fmhL = 0.5*minmod_dx(gp, gmhL, gmhL2);
                    double fmhR = -0.5*minmod_dx(gphR, gp, gmhL);
                    qiphL[alpha] = gphL + fphL;
                    qiphR[alpha] = gphR + fphR;
                    qimhL[alpha] = gmhL + fmhL;
                    qimhR[alpha] = gmhR + fmhR;
                }
                // for each direction, reconstruct half-way cells
                // reconstruct e, rhob, and u[4] for half way cells
                flag = ReconstIt_velocity_Newton(
                                grid_array_hL, tau, qiphL, grid_array[idx]);
                aiphL = MaxSpeed(tau, direc, grid_array_hL);

                flag *= ReconstIt_velocity_Newton(
                                grid_array_hR, tau, qiphR, grid_array[idx]);
                aiphR = MaxSpeed(tau, direc, grid_array_hR);
                aiph = maxi(aiphL, aiphR);
                for (int alpha = 0; alpha < 5; alpha++) {
                    double FiphL = tau_fac*get_TJb_new(grid_array_hL,
                                                       alpha, direc);
                    double FiphR = tau_fac*get_TJb_new(grid_array_hR,
                                                        alpha, direc);
                    // KT: H_{j+1/2} = (f(u^+_{j+1/2}) + f(u^-_{j+1/2})/2
                    //              - a_{j+1/2}(u_{j+1/2}^+ - u^-_{j+1/2})/2
                    double Fiph = 0.5*((FiphL + FiphR)
                                        - aiph*(qiphR[alpha] - qiphL[alpha]));

                    rhs[alpha] -= Fiph/0.1*0.01;
                }

                flag *= ReconstIt_velocity_Newton(grid_array_hL, tau,
                                                     qimhL, grid_array[idx]);
                aimhL = MaxSpeed(tau, direc, grid_array_hL);

                flag *= ReconstIt_velocity_Newton(grid_array_hR, tau,
                                                     qimhR, grid_array[idx]);
                aimhR = MaxSpeed(tau, direc, grid_array_hR);
                aimh = maxi(aimhL, aimhR);

                for (int alpha = 0; alpha < 5; alpha++) {
                    double FimhL = tau_fac*get_TJb_new(grid_array_hL,
                                                       alpha, direc);
                    double FimhR = tau_fac*get_TJb_new(grid_array_hR,
                                                       alpha, direc);
                    // KT: H_{j+1/2} = (f(u^+_{j+1/2}) + f(u^-_{j+1/2})/2
                    //              - a_{j+1/2}(u_{j+1/2}^+ - u^-_{j+1/2})/2
                    double Fimh = 0.5*((FimhL + FimhR)
                                        - aimh*(qimhR[alpha] - qimhL[alpha]));
                    rhs[alpha] += Fimh/0.1*0.01;
                }
                //cout << "y-direction" << endl;
                
                // eta-direction
                direc = 3;
                tau_fac = 1.0;
                for (int alpha = 0; alpha < 5; alpha++) {
                    double gp = qi_array[idx][alpha];
                    double gphL = qi_array[idx][alpha];
                    double gmhR = qi_array[idx][alpha];

                    double gphR, gmhL, gphR2, gmhL2;
                    if (k + 1 < n_cell_eta) {
                        int idx_p_1 = j + i*n_cell_y + (k+1)*n_cell_x*n_cell_y;
                        gphR = qi_array[idx_p_1][alpha];
                    } else {
                        int idx_p_1 = 4*j + 4*i*n_cell_y + 2;
                        gphR = qi_nbr_eta[idx_p_1][alpha];
                    }
                    if (k - 1 >= 0) {
                        int idx_m_1 = j + i*n_cell_y + (k-1)*n_cell_x*n_cell_y;
                        gmhL = qi_array[idx_m_1][alpha];
                    } else {
                        int idx_m_1 = 4*j + 4*i*n_cell_y + 1;
                        gmhL = qi_nbr_eta[idx_m_1][alpha];
                    }
                    if (k + 2 < n_cell_eta) {
                        int idx_p_2 = j + i*n_cell_y + (k+2)*n_cell_x*n_cell_y;
                        gphR2 = qi_array[idx_p_2][alpha];
                    } else {
                        int idx_p_2 = 4*j + 4*i*n_cell_y + 4 + k - n_cell_eta;
                        gphR2 = qi_nbr_eta[idx_p_2][alpha];
                    }
                    if (k - 2 >= 0) {
                        int idx_m_2 = j + i*n_cell_y + (k-2)*n_cell_x*n_cell_y;
                        gmhL2 = qi_array[idx_m_2][alpha];
                    } else {
                        int idx_m_2 = 4*j + 4*i*n_cell_y + k;
                        gmhL2 = qi_nbr_eta[idx_m_2][alpha];
                    }

                    double fphL = 0.5*minmod_dx(gphR, gp, gmhL);
                    double fphR = -0.5*minmod_dx(gphR2, gphR, gp);
                    double fmhL = 0.5*minmod_dx(gp, gmhL, gmhL2);
                    double fmhR = -0.5*minmod_dx(gphR, gp, gmhL);
                    qiphL[alpha] = gphL + fphL;
                    qiphR[alpha] = gphR + fphR;
                    qimhL[alpha] = gmhL + fmhL;
                    qimhR[alpha] = gmhR + fmhR;
                }
                // for each direction, reconstruct half-way cells
                // reconstruct e, rhob, and u[4] for half way cells
                flag = ReconstIt_velocity_Newton(
                                grid_array_hL, tau, qiphL, grid_array[idx]);
                aiphL = MaxSpeed(tau, direc, grid_array_hL);

                flag *= ReconstIt_velocity_Newton(
                                grid_array_hR, tau, qiphR, grid_array[idx]);
                aiphR = MaxSpeed(tau, direc, grid_array_hR);
                aiph = maxi(aiphL, aiphR);
                for (int alpha = 0; alpha < 5; alpha++) {
                    double FiphL = tau_fac*get_TJb_new(grid_array_hL,
                                                       alpha, direc);
                    double FiphR = tau_fac*get_TJb_new(grid_array_hR,
                                                        alpha, direc);
                    // KT: H_{j+1/2} = (f(u^+_{j+1/2}) + f(u^-_{j+1/2})/2
                    //              - a_{j+1/2}(u_{j+1/2}^+ - u^-_{j+1/2})/2
                    double Fiph = 0.5*((FiphL + FiphR)
                                        - aiph*(qiphR[alpha] - qiphL[alpha]));

                    rhs[alpha] -= Fiph/0.1*0.01;
                }

                flag *= ReconstIt_velocity_Newton(grid_array_hL, tau,
                                                     qimhL, grid_array[idx]);
                aimhL = MaxSpeed(tau, direc, grid_array_hL);

                flag *= ReconstIt_velocity_Newton(grid_array_hR, tau,
                                                     qimhR, grid_array[idx]);
                aimhR = MaxSpeed(tau, direc, grid_array_hR);
                aimh = maxi(aimhL, aimhR);

                for (int alpha = 0; alpha < 5; alpha++) {
                    double FimhL = tau_fac*get_TJb_new(grid_array_hL,
                                                       alpha, direc);
                    double FimhR = tau_fac*get_TJb_new(grid_array_hR,
                                                       alpha, direc);
                    // KT: H_{j+1/2} = (f(u^+_{j+1/2}) + f(u^-_{j+1/2})/2
                    //              - a_{j+1/2}(u_{j+1/2}^+ - u^-_{j+1/2})/2
                    double Fimh = 0.5*((FimhL + FimhR)
                                        - aimh*(qimhR[alpha] - qimhL[alpha]));
                    rhs[alpha] += Fimh/0.1*0.01;
                }
                //cout << "eta-direction" << endl;

                // geometric terms
                rhs[0] -= (get_TJb_new(grid_array[idx], 3, 3)
                           *0.01);
                rhs[3] -= (get_TJb_new(grid_array[idx], 3, 0)
                           *0.01);
                
                for (int alpha = 0; alpha < 5; alpha++) {
                    qi_array_new[idx][alpha] = qi_array[idx][alpha] + rhs[alpha];
                }
            }
        }
    }

}

int test::ReconstIt_velocity_Newton(
    double *grid_array, double tau, double *uq, double *grid_array_p) {
    double abs_err = 1e-9;
    double rel_err = 1e-8;
    int max_iter = 100;
    double v_critical = 0.53;
    /* prepare for the iteration */
    /* uq = qiphL, qiphR, etc 
       qiphL[alpha] means, for instance, TJ[alpha][0] 
       in the cell at x+dx/2 calculated from the left */
    
    /* uq are the conserved charges. That is, the ones appearing in
       d_tau (Ttautau/tau) + d_eta(Ttaueta/tau) + d_perp(Tperptau) = -Tetaeta
       d_tau (Ttaueta) + d_eta(Tetaeta) + d_v(tau Tveta) = -Ttaueta/tau 
       d_tau (Ttauv) + d_eta(Tetav) + d_w(tau Twv) = 0
       d_tau (Jtau) + d_eta Jeta + d_perp(tau Jperp) = 0 */
    
    /* q[0] = Ttautau/tau, q[1] = Ttaux, q[2] = Ttauy, q[3] = Ttaueta
       q[4] = Jtau */
    /* uq = qiphL, qiphR, qimhL, qimhR, qirk */

    double K00 = (uq[1]*uq[1] + uq[2]*uq[2] + uq[3]*uq[3])/(tau*tau);
    double M = sqrt(K00);
    double T00 = uq[0]/tau;
    double J0 = uq[4]/tau;

    if ((T00 < abs_err) || ((T00 - K00/T00) < 0.0)) {
        revert_grid(grid_array, grid_array_p);
    }/* if t00-k00/t00 < 0.0 */

    double u0, u1, u2, u3, epsilon, pressure, rhob;

    double u0_guess = 1./sqrt(1. - grid_array_p[1]*grid_array_p[1]
                              - grid_array_p[2]*grid_array_p[2]
                              - grid_array_p[3]*grid_array_p[3]);
    double v_guess = sqrt(1. - 1./(u0_guess*u0_guess + 1e-15));
    //if (isnan(v_guess)) {
    //    v_guess = 0.0;
    //}
    int v_status = 1;
    int iter = 0;
    double rel_error_v = 10.0;
    double v_next = 1.0;
    double v_prev = v_guess;
    double abs_error_v = reconst_velocity_f_Newton(v_prev, T00, M, J0);
    do {
        iter++;
        v_next = (v_prev
                  - (abs_error_v/reconst_velocity_df(v_prev, T00, M, J0)));
        if (v_next < 0.0) {
            v_next = 0.0 + 1e-10;
        } else if (v_next > 1.0) {
            v_next = 1.0 - 1e-10;
        }
        abs_error_v = reconst_velocity_f_Newton(v_next, T00, M, J0);
        rel_error_v = 2.*abs_error_v/(v_next + v_prev + 1e-15);
        v_prev = v_next;
        if (iter > max_iter) {
            v_status = 0;
            break;
        }
    } while (abs_error_v > abs_err && rel_error_v > rel_err);

    double v_solution;
    if (v_status == 1) {
        v_solution = v_next;
    } else {
        revert_grid(grid_array, grid_array_p);
    }/* if iteration is unsuccessful, revert */
   
    // for large velocity, solve u0
    double u0_solution = 1.0;
    if (v_solution > v_critical) {
        double u0_prev = 1./sqrt(1. - v_solution*v_solution);
        int u0_status = 1;
        iter = 0;
        double u0_next;
        double abs_error_u0 = reconst_u0_f_Newton(u0_prev, T00, K00, M, J0);
        double rel_error_u0 = 1.0;
        do {
            iter++;
            u0_next = (u0_prev
                       - abs_error_u0/reconst_u0_df(u0_prev, T00, K00, M, J0));
            abs_error_u0 = reconst_u0_f_Newton(u0_next, T00, K00, M, J0);
            rel_error_u0 = 2.*abs_error_u0/(u0_next + u0_prev + 1e-15);
            u0_prev = u0_next;
            if (iter > max_iter) {
                u0_status = 0;
                break;
            }
        } while (abs_error_u0 > abs_err && rel_error_u0 > rel_err);

        if (u0_status == 1) {
            u0_solution = u0_next;
        } else {
        }  // if iteration is unsuccessful, revert
    }

    // successfully found velocity, now update everything else
    if (v_solution < v_critical) {
        u0 = 1./(sqrt(1. - v_solution*v_solution) + v_solution*abs_err);
        epsilon = T00 - v_solution*sqrt(K00);
        rhob = J0/u0;
    } else {
        u0 = u0_solution;
        epsilon = T00 - sqrt((1. - 1./(u0_solution*u0_solution))*K00);
        rhob = J0/u0_solution;
    }

    double check_u0_var = (fabs(u0 - u0_guess)/u0_guess);
    if (check_u0_var > 100.) {
        revert_grid(grid_array, grid_array_p);
    }

    grid_array[0] = epsilon;
    grid_array[4] = rhob;

    pressure = get_pressure(epsilon, rhob);

    // individual components of velocity
    double velocity_inverse_factor = u0/(T00 + pressure);

    double u_max = 242582597.70489514; // cosh(20)
    //remove if for speed
    if(u0 > u_max) {
        revert_grid(grid_array, grid_array_p);
    } else {
        u1 = uq[1]*velocity_inverse_factor/tau; 
        u2 = uq[2]*velocity_inverse_factor/tau; 
        u3 = uq[3]*velocity_inverse_factor/tau;
    }

    // Correcting normalization of 4-velocity
    double temp_usq = u0*u0 - u1*u1 - u2*u2 - u3*u3;
    // Correct velocity when unitarity is not satisfied to numerical accuracy
    if (fabs(temp_usq - 1.0) > abs_err) {
        // If the deviation is too large, exit MUSIC
        if (fabs(temp_usq - 1.0) > 0.1*u0) {
            revert_grid(grid_array, grid_array_p);
        }
        // Rescaling spatial components of velocity so that unitarity 
        // is exactly satisfied (u[0] is not modified)
        double scalef = sqrt((u0*u0 - 1.0)
                             /(u1*u1 + u2*u2 + u3*u3 + abs_err));
        u1 *= scalef;
        u2 *= scalef;
        u3 *= scalef;
    }  // if u^mu u_\mu != 1 
    // End: Correcting normalization of 4-velocity
   
    grid_array[1] = u1/u0;
    grid_array[2] = u2/u0;
    grid_array[3] = u3/u0;

    return(1);  /* on successful execution */
}

void test::revert_grid(double *grid_array, double *grid_prev) {
    for (int i = 0; i < 5; i++) {
        grid_array[i] = grid_prev[i];
    }
}

double test::get_pressure(double e_local, double rhob) {
    double p = e_local/3.;
    return(p);
}

double test::MaxSpeed(double tau, int direc, double *grid_array) {
    //grid_p = grid_p_h_L, grid_p_h_R, grid_m_h_L, grid_m_h_R
    //these are reconstructed by Reconst, grid_array = [e, v^i, rhob]
    //double utau = (grid_p->u[0][0]);
    double utau = 1./sqrt(1. - grid_array[1]*grid_array[1]
                             - grid_array[2]*grid_array[2]
                             - grid_array[3]*grid_array[3]);
    double utau2 = utau*utau;
    //double ux = fabs((grid_p->u[0][direc]));
    double ux = fabs(utau*grid_array[direc]);
    double ux2 = ux*ux;
    double ut2mux2 = utau2 - ux2;
  
    //double eps = grid_p->epsilon;
    //double rhob = grid_p->rhob;
    double eps = grid_array[0];
    double rhob = grid_array[4];
  
    double vs2 = get_cs2(eps, rhob);

    double den = utau2*(1. - vs2) + vs2;
    double num_temp_sqrt = (ut2mux2 - (ut2mux2 - 1.)*vs2)*vs2;
    double num;
    if (num_temp_sqrt >= 0) {
        num = utau*ux*(1. - vs2) + sqrt(num_temp_sqrt);
    } else {
        double dpde = p_e_func(eps, rhob);
        double p = get_pressure(eps, rhob);
        double h = p + eps;
        if (dpde < 0.001) {
            num = (sqrt(-(h*dpde*h*(dpde*(-1.0 + ut2mux2) - ut2mux2))) 
                   - h*(-1.0 + dpde)*utau*ux);
        } else {
            num = 0.0;
        }
    }
    
    double f = num/(den + 1e-15);
    if (f > 1.0) {
        f = 1.0;
    }
    if (f < ux/utau) {
        f = ux/utau;
    }

    if (direc == 3) {
        f /= tau;
    }

    return(f);
}

double test::minmod_dx(double up1, double u, double um1) {
    double theta_flux = 1.8;
    double diffup = (up1 - u)*theta_flux;
    double diffdown = (u - um1)*theta_flux;
    double diffmid = (up1 - um1)*0.5;

    double tempf;
    if ( (diffup > 0.0) && (diffdown > 0.0) && (diffmid > 0.0) ) {
        tempf = mini(diffdown, diffmid);
        return mini(diffup, tempf);
    } else if ( (diffup < 0.0) && (diffdown < 0.0) && (diffmid < 0.0) ) {
        tempf = maxi(diffdown, diffmid);
        return maxi(diffup, tempf);
    } else {
      return 0.0;
    }
}/* minmod_dx */

double test::get_TJb_new(double *grid_array, int mu, int nu) {
    double rhob = grid_array[4];
    double gamma = 1./sqrt(1. - grid_array[1]*grid_array[1]
                              - grid_array[2]*grid_array[2]
                              - grid_array[3]*grid_array[3]);
    double u_nu = gamma;
    if (nu > 0) {
        u_nu *= grid_array[nu];
    }

    if (mu == 4) {
        double J_nu = rhob*u_nu;
        return(J_nu);
    } else if (mu < 4) {
        double e = grid_array[0];
        double gfac = 0.0;
        double u_mu = gamma;
        if (mu == nu) {
            u_mu = u_nu;
            if (mu == 0) {
                gfac = -1.0;
            } else {
                gfac = 1.0;
            }
        } else {
            if (mu > 0) {
                u_mu *= grid_array[mu];
            }
        }
        double pressure = get_pressure(e, rhob);
        double T_munu = (e + pressure)*u_mu*u_nu + pressure*gfac;
        return(T_munu);
    } else {
        return(0.0);
    }
}

double test::reconst_velocity_f(double v, double T00, double M,
                                   double J0) {
    // this function returns f(v) = M/(M0 + P)
    double epsilon = T00 - v*M;
    double rho = J0*sqrt(1 - v*v);
   
    double pressure = get_pressure(epsilon, rho);
    double fv = M/(T00 + pressure);
    return(fv);
}

double test::reconst_velocity_f_Newton(double v, double T00, double M,
                                          double J0) {
    double fv = v - reconst_velocity_f(v, T00, M, J0);
    return(fv);
}

double test::reconst_velocity_df(double v, double T00, double M,
                                    double J0) {
    // this function returns df'(v)/dv where f' = v - f(v)
    double epsilon = T00 - v*M;
    double temp = sqrt(1. - v*v);
    double rho = J0*temp;
    double temp2 = v/temp;
   
    double pressure = get_pressure(epsilon, rho);
    double dPde = p_e_func(epsilon, rho);
    double dPdrho = p_rho_func(epsilon, rho);
    
    double temp1 = T00 + pressure;

    double dfdv = 1. - M/(temp1*temp1)*(M*dPde + J0*temp2*dPdrho);
    return(dfdv);
}

double test::reconst_u0_f(double u0, double T00, double K00, double M,
                             double J0) {
    // this function returns f(u0) = (M0+P)/sqrt((M0+P)^2 - M^2)
    double epsilon = T00 - sqrt(1. - 1./u0/u0)*M;
    double rho = J0/u0;
    
    double pressure = get_pressure(epsilon, rho);
    double fu = (T00 + pressure)/sqrt((T00 + pressure)*(T00 + pressure) - K00);
    return(fu);
}

double test::reconst_u0_f_Newton(double u0, double T00, double K00,
                                    double M, double J0) {
    // this function returns f(u0) = u0 - (M0+P)/sqrt((M0+P)^2 - M^2)
    double fu = u0 - reconst_u0_f(u0, T00, K00, M, J0);
    return(fu);
}

double test::reconst_u0_df(double u0, double T00, double K00, double M,
                              double J0) {
    // this function returns df'/du0 where f'(u0) = u0 - f(u0)
    double v = sqrt(1. - 1./(u0*u0));
    double epsilon = T00 - v*M;
    double rho = J0/u0;
    double dedu0 = - M/(u0*u0*u0*v);
    double drhodu0 = - J0/(u0*u0);
    
    double pressure = get_pressure(epsilon, rho);
    double dPde = p_e_func(epsilon, rho);
    double dPdrho = p_rho_func(epsilon, rho);

    double denorm = pow(((T00 + pressure)*(T00 + pressure) - K00), 1.5);
    double dfdu0 = 1. + (dedu0*dPde + drhodu0*dPdrho)*K00/denorm;
    return(dfdu0);
}

double test::get_cs2(double e_local, double rhob) {
    double cs2 = 1./3.;
    return(cs2);
}

double test::p_e_func(double e_local, double rhob) {
    double dPde = 1./3.;
    return(dPde);
}

double test::p_rho_func(double e_local, double rhob) {
    double dPdrho = 0.0;
    return(dPdrho);
}

double test::get_mu(double e_local, double rhob) {
    double mu = 0.0;
    return(mu);
}

void test::update_grid_cell(double grid_array[][5], Field *hydro_fields, int rk_flag,
                               int ieta, int ix, int iy,
                               int n_cell_eta, int n_cell_x, int n_cell_y) {
    for (int k = 0; k < n_cell_eta; k++) {
        int idx_ieta = get_min(ieta + k, 0);
        for (int i = 0; i < n_cell_x; i++) {
            int idx_ix = get_min(ix + i, 200);
            for (int j = 0; j < n_cell_y; j++) {
                int idx_iy = get_min(iy + j, 200);
                int field_idx = (idx_iy + idx_ix*(201)
                                 + idx_ieta*(201)*(201));
                int idx = j + i*n_cell_y + k*n_cell_x*n_cell_y;
                update_grid_array_to_hydro_fields(
                        grid_array[idx], hydro_fields, field_idx, rk_flag);
            }
        }
    }
}           

void test::update_grid_array_to_hydro_fields(
            double *grid_array, Field *hydro_fields, int idx, int rk_flag) {
    double gamma = 1./sqrt(1. - grid_array[1]*grid_array[1]
                              - grid_array[2]*grid_array[2]
                              - grid_array[3]*grid_array[3]);
    if (rk_flag == 0) {
        hydro_fields->e_rk1[idx] = grid_array[0];
        hydro_fields->rhob_rk1[idx] = grid_array[4];
        hydro_fields->u_rk1[idx][0] = gamma;
        for (int i = 1; i < 4; i++) {
            hydro_fields->u_rk1[idx][i] = gamma*grid_array[i];
        }
    } else {
        hydro_fields->e_rk0[idx] = grid_array[0];
        hydro_fields->rhob_rk0[idx] = grid_array[4];
        hydro_fields->u_rk0[idx][0] = gamma;
        for (int i = 1; i < 4; i++) {
            hydro_fields->u_rk0[idx][i] = gamma*grid_array[i];
        }
    }
}
