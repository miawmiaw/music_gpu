// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include "./util.h"
#include "./grid.h"
#include "./data.h"
#include "./eos.h"
#include "./dissipative.h"
#include<vector>

using namespace std;

Diss::Diss(EOS *eosIn, InitData* DATA_in) {
    eos = eosIn;
    DATA_ptr = DATA_in;
    minmod = new Minmod(DATA_in);
    util = new Util;
}

// destructor
Diss::~Diss() {
    delete minmod;
    delete util;
}


void Diss::MakeWSource(double tau, double **qi_array,
                       int n_cell_eta, int n_cell_x, int n_cell_y,
                       double **vis_array,
                       double **vis_nbr_tau, double **vis_nbr_x,
                       double **vis_nbr_y, double **vis_nbr_eta,
                       double **qi_array_new) {
//! calculate d_m (tau W^{m,alpha}) + (geom source terms)
//! partial_tau W^tau alpha
//! this is partial_tau evaluated at tau
//! this is the first step. so rk_flag = 0
//! change: alpha first which is the case
//!         for everywhere else. also, this change is necessary
//!         to use Wmunu[rk_flag][4][mu] as the dissipative baryon current

    double shear_on, bulk_on;
    if (DATA_ptr->turn_on_shear)
        shear_on = 1.0;
    else
        shear_on = 0.0;

    if (DATA_ptr->turn_on_bulk)
        bulk_on = 1.0;
    else
        bulk_on = 0.0;

    int alpha_max = 5;
    if (DATA_ptr->turn_on_diff == 0) {
        alpha_max = 4;
    }
    for (int k = 0; k < n_cell_eta; k++) {
        for (int i = 0; i < n_cell_x; i++) {
            for (int j = 0; j < n_cell_y; j++) {
                int idx = j + i*n_cell_y + k*n_cell_x*n_cell_y;
                for (int alpha = 0; alpha < alpha_max; alpha++) {
                    // dW/dtau
                    // backward time derivative (first order is more stable)
                    int idx_1d_alpha0 = util->map_2d_idx_to_1d(alpha, 0);
                    double dWdtau;
                    dWdtau = ((vis_array[idx][idx_1d_alpha0]
                               - vis_nbr_tau[idx][idx_1d_alpha0])
                              /DATA_ptr->delta_tau);

                    // bulk pressure term
                    double dPidtau = 0.0;
                    double Pi_alpha0 = 0.0;
                    if (alpha < 4 && DATA_ptr->turn_on_bulk == 1) {
                        double gfac = (alpha == 0 ? -1.0 : 0.0);
                        Pi_alpha0 = (vis_array[idx][14]
                                     *(gfac + vis_array[idx][15+alpha]
                                              *vis_array[idx][15]));

                        dPidtau = (Pi_alpha0
                                   - vis_nbr_tau[idx][14]
                                     *(gfac + vis_nbr_tau[idx][alpha+15]
                                              *vis_nbr_tau[idx][15]));
                    }

                    // use central difference to preserve
                    // the conservation law exactly
                    int idx_1d;
                    int idx_p_1, idx_m_1;
                    double dWdx_perp = 0.0;
                    double dPidx_perp = 0.0;

                    double sgp1, sgm1, bgp1, bgm1;
                    // x-direction
                    idx_1d = util->map_2d_idx_to_1d(alpha, 1);
                    if (i + 1 < n_cell_x) {
                        idx_p_1 = j + (i+1)*n_cell_y + k*n_cell_x*n_cell_y;
                        sgp1 = vis_array[idx_p_1][idx_1d];
                    } else {
                        idx_p_1 = 4*j + k*4*n_cell_y + 2;
                        sgp1 = vis_nbr_x[idx_p_1][idx_1d];
                    }
                    if (i - 1 >= 0) {
                        idx_m_1 = j + (i-1)*n_cell_y + k*n_cell_x*n_cell_y;
                        sgm1 = vis_array[idx_m_1][idx_1d];
                    } else {
                        idx_m_1 = 4*j + k*4*n_cell_y + 1;
                        sgm1 = vis_nbr_x[idx_m_1][idx_1d];
                    }
                    dWdx_perp += (sgp1 - sgm1)/(2.*DATA_ptr->delta_x);
                    if (alpha < 4 && DATA_ptr->turn_on_bulk == 1) {
                        double gfac1 = (alpha == 1 ? 1.0 : 0.0);
                        if (i + 1 < n_cell_x) {
                            bgp1 = (vis_array[idx_p_1][14]
                                        *(gfac1 + vis_array[idx_p_1][15+alpha]
                                                  *vis_nbr_x[idx_p_1][16]));
                        } else {
                            bgp1 = (vis_nbr_x[idx_p_1][14]
                                        *(gfac1 + vis_nbr_x[idx_p_1][15+alpha]
                                                  *vis_nbr_x[idx_p_1][16]));
                        }
                        if (i - 1 >= 0) {
                            bgm1 = (vis_array[idx_m_1][14]
                                        *(gfac1 + vis_array[idx_m_1][15+alpha]
                                                  *vis_array[idx_m_1][16]));
                        } else {
                            bgm1 = (vis_nbr_x[idx_m_1][14]
                                        *(gfac1 + vis_nbr_x[idx_m_1][15+alpha]
                                                  *vis_nbr_x[idx_m_1][16]));
                        }
                        dPidx_perp += (bgp1 - bgm1)/(2.*DATA_ptr->delta_x);
                    }
                    // y-direction
                    idx_1d = util->map_2d_idx_to_1d(alpha, 2);
                    if (j + 1 < n_cell_y) {
                        idx_p_1 = j + 1 + i*n_cell_y + k*n_cell_x*n_cell_y;
                        sgp1 = vis_array[idx_p_1][idx_1d];
                    } else {
                        idx_p_1 = 4*i + 4*k*n_cell_x + 2;
                        sgp1 = vis_nbr_y[idx_p_1][idx_1d];
                    }
                    if (j - 1 >= 0) {
                        idx_m_1 = j - 1 + i*n_cell_y + k*n_cell_x*n_cell_y;
                        sgm1 = vis_array[idx_m_1][idx_1d];
                    } else {
                        idx_m_1 = 4*i + 4*k*n_cell_x + 1;
                        sgm1 = vis_nbr_y[idx_m_1][idx_1d];
                    }
                    dWdx_perp += (sgp1 - sgm1)/(2.*DATA_ptr->delta_x);
                    if (alpha < 4 && DATA_ptr->turn_on_bulk == 1) {
                        double gfac1 = (alpha == 2 ? 1.0 : 0.0);
                        if (j + 1 < n_cell_x) {
                            bgp1 = (vis_array[idx_p_1][14]
                                        *(gfac1 + vis_nbr_y[idx_p_1][15+alpha]
                                                  *vis_nbr_y[idx_p_1][17]));
                        } else {
                            bgp1 = (vis_nbr_y[idx_p_1][14]
                                        *(gfac1 + vis_nbr_y[idx_p_1][15+alpha]
                                                  *vis_nbr_y[idx_p_1][17]));
                        }
                        if (j - 1 >= 0) {
                            bgm1 = (vis_array[idx_m_1][14]
                                        *(gfac1 + vis_array[idx_m_1][15+alpha]
                                                  *vis_array[idx_m_1][17]));
                        } else {
                            bgm1 = (vis_nbr_y[idx_m_1][14]
                                        *(gfac1 + vis_nbr_y[idx_m_1][15+alpha]
                                                  *vis_nbr_y[idx_m_1][17]));
                        }
                        dPidx_perp += (bgp1 - bgm1)/(2.*DATA_ptr->delta_x);
                    }

                    // eta-direction
                    double taufactor = tau;
                    double dWdeta = 0.0;
                    double dPideta = 0.0;
                    idx_1d = util->map_2d_idx_to_1d(alpha, 3);
                    if (k + 1 < n_cell_eta) {
                        idx_p_1 = j + i*n_cell_y + (k+1)*n_cell_x*n_cell_y;
                        sgp1 = vis_array[idx_p_1][idx_1d];
                    } else {
                        idx_p_1 = 4*j + 4*i*n_cell_y + 2;
                        sgp1 = vis_nbr_eta[idx_p_1][idx_1d];
                    }
                    if (k - 1 >= 0) {
                        idx_m_1 = j + i*n_cell_y + (k-1)*n_cell_x*n_cell_y;
                        sgm1 = vis_array[idx_m_1][idx_1d];
                    } else {
                        idx_m_1 = 4*j + 4*i*n_cell_y + 1;
                        sgm1 = vis_nbr_eta[idx_m_1][idx_1d];
                    }
                    dWdeta = (sgp1 - sgm1)/(2.*DATA_ptr->delta_eta*taufactor);
                    if (alpha < 4 && DATA_ptr->turn_on_bulk == 1) {
                        double gfac3 = (alpha == 3 ? 1.0 : 0.0);
                        if (k + 1 < n_cell_eta) {
                            bgp1 = (vis_array[idx_p_1][14]
                                       *(gfac3 + vis_array[idx_p_1][15+alpha]
                                                 *vis_array[idx_p_1][18]));
                        } else {
                            bgp1 = (vis_nbr_eta[idx_p_1][14]
                                       *(gfac3 + vis_nbr_eta[idx_p_1][15+alpha]
                                                 *vis_nbr_eta[idx_p_1][18]));
                        }
                        if (k - 1 >= 0) {
                            bgm1 = (vis_array[idx_m_1][14]
                                       *(gfac3 + vis_array[idx_m_1][15+alpha]
                                                 *vis_array[idx_m_1][18]));
                        } else {
                            bgm1 = (vis_nbr_eta[idx_m_1][14]
                                       *(gfac3 + vis_nbr_eta[idx_m_1][15+alpha]
                                                 *vis_nbr_eta[idx_m_1][18]));
                        }
                        dPideta = ((bgp1 - bgm1)
                                   /(2.*DATA_ptr->delta_eta*taufactor));
                    }

                    // partial_m (tau W^mn) = W^0n + tau partial_m W^mn
                    double sf = (tau*(dWdtau + dWdx_perp + dWdeta)
                                 + vis_array[idx][idx_1d_alpha0]);
                    double bf = (tau*(dPidtau + dPidx_perp + dPideta)
                                 + Pi_alpha0);

                    // sources due to coordinate transform
                    // this is added to partial_m W^mn
                    if (alpha == 0) {
                        sf += vis_array[idx][9];
                        bf += vis_array[idx][14]*(1.0 + vis_array[idx][18]
                                                        *vis_array[idx][18]);
                    }
                    if (alpha == 3) {
                        sf += vis_array[idx][3];
                        bf += vis_array[idx][14]*(vis_array[idx][15]
                                                  *vis_array[idx][18]);
                    }

                    double result = 0.0;
                    if (alpha < 4) {
                        result = (sf*shear_on + bf*bulk_on);
                    } else if (alpha == 4) {
                        result = sf;
                    }
                    qi_array_new[idx][alpha] -= result*(DATA_ptr->delta_tau);
                }
            }
        }
    }
}


double Diss::Make_uWSource(double tau, int n_cell_eta, int n_cell_x,
                           int n_cell_y,
                           double **vis_array, double **velocity_array,
                           double **grid_array, double **vis_array_new) {
    if (DATA_ptr->turn_on_shear == 0) {
        return(0.0);
    }

    int include_WWterm = 1;
    int include_Vorticity_term = 0;
    int include_Wsigma_term = 1;
    if (DATA_ptr->Initial_profile == 0) {
        include_WWterm = 0;
        include_Wsigma_term = 0;
        include_Vorticity_term = 0;
    }

    double sigma[4][4];
    double Wmunu[4][4];
    for (int k = 0; k < n_cell_eta; k++) {
        for (int i = 0; i < n_cell_x; i++) {
            for (int j = 0; j < n_cell_y; j++) {
                int idx = j + i*n_cell_y + k*n_cell_x*n_cell_y;

                for (int a = 0; a < 4; a++) {
                    for (int b = a; b < 4; b++) {
                        int idx_1d = util->map_2d_idx_to_1d(a, b);
                        //Wmunu[a][b] = grid_pt->Wmunu[rk_flag][idx_1d];
                        //sigma[a][b] = sigma_1d[idx_1d];
                        Wmunu[a][b] = vis_array[idx][idx_1d];
                        sigma[a][b] = velocity_array[idx][6+idx_1d];
                    }
                }
                for (int a = 0; a < 4; a++) {
                    for (int b = a+1; b < 4; b++) {
                        Wmunu[b][a] = Wmunu[a][b];
                        sigma[b][a] = sigma[a][b];
                    }
                }

                // Useful variables to define
                double epsilon = grid_array[idx][0];
                double rhob = grid_array[idx][4];

                double T = eos->get_temperature(epsilon, rhob);

                double shear_to_s = 0.0;
                if (DATA_ptr->T_dependent_shear_to_s == 1) {
                    shear_to_s = get_temperature_dependent_eta_s(T);
                } else {
                    shear_to_s = DATA_ptr->shear_to_s;
                }


                //  Defining transport coefficients  
                double pressure = eos->get_pressure(epsilon, rhob);
                double shear = (shear_to_s)*(epsilon + pressure)/(T + 1e-15);
                double tau_pi = 5.0*shear/(epsilon + pressure + 1e-15);
                if (tau_pi < DATA_ptr->delta_tau) {
                    tau_pi = DATA_ptr->delta_tau;
                }

                // transport coefficient for nonlinear terms
                // -- shear only terms -- 4Mar2013
                // transport coefficients of a massless gas of
                // single component particles
                double transport_coefficient  = 9./70.*tau_pi/shear*(4./5.);
                double transport_coefficient2 = 4./3.*tau_pi;
                double transport_coefficient3 = 10./7.*tau_pi;
                double transport_coefficient4 = 2.*tau_pi;

                // transport coefficient for nonlinear terms
                // -- coupling to bulk viscous pressure -- 4Mar2013
                // transport coefficients not yet known -- fixed to zero
                double transport_coefficient_b  = 6./5.*tau_pi;
                double transport_coefficient2_b = 0.;

                // This source has many terms
                // everything in the 1/(tau_pi) piece is here
                // third step in the split-operator time evol
                //  use Wmunu[rk_flag] and u[rk_flag] with rk_flag = 0

                // Wmunu + transport_coefficient2*Wmunu*theta

                for (int mu = 1; mu < 4; mu++) {
                    for (int nu = mu; nu < 4; nu++) {
                        int idx_1d = util->map_2d_idx_to_1d(mu, nu);
                        // full term is
                        //- (1.0 + transport_coefficient2*theta_local)
                        double tempf = (
                        - (1.0 + transport_coefficient2*velocity_array[idx][0])
                          *(Wmunu[mu][nu]));

                        // Navier-Stokes Term -- -2.*shear*sigma^munu
                        // full Navier-Stokes term is
                        // sign changes according to metric sign convention
                        double NS_term = - 2.*shear*sigma[mu][nu];

                        // Vorticity Term
                        double Vorticity_term = 0.0;
                        // for future
                        // remember: dUsup[m][n] = partial^n u^m  ///
                        // remember:  a[n]  =  u^m*partial_m u^n  ///
                        //if (include_Vorticity_term == 1) {
                        //    double term1_Vorticity;
                        //    double omega[4][4];
                        //    for (a = 0; a < 4; a++) {
                        //        for (b = 0; b <4; b++) {
                        //            omega[a][b] = (
                        //                (grid_pt->dUsup[0][a][b]
                        //                 - grid_pt->dUsup[0][b][a])/2.
                        //                + ueta/tau/2.*(DATA->gmunu[a][0]*DATA->gmunu[b][3]
                        //                               - DATA->gmunu[b][0]*DATA->gmunu[a][3])
                        //                - ueta*gamma/tau/2.
                        //                  *(DATA->gmunu[a][3]*grid_pt->u[rk_flag][b]
                        //                    - DATA->gmunu[b][3]*grid_pt->u[rk_flag][a])
                        //                + ueta*ueta/tau/2.
                        //                  *(DATA->gmunu[a][0]*grid_pt->u[rk_flag][b]
                        //                     - DATA->gmunu[b][0]*grid_pt->u[rk_flag][a])
                        //                + (grid_pt->u[rk_flag][a]*a_local[b]
                        //                   - grid_pt->u[rk_flag][b]*a_local[a])/2.);
                        //        }
                        //    }
                        //    term1_Vorticity = (- Wmunu[mu][0]*omega[nu][0]
                        //                       - Wmunu[nu][0]*omega[mu][0]
                        //                       + Wmunu[mu][1]*omega[nu][1]
                        //                       + Wmunu[nu][1]*omega[mu][1]
                        //                       + Wmunu[mu][2]*omega[nu][2]
                        //                       + Wmunu[nu][2]*omega[mu][2]
                        //                       + Wmunu[mu][3]*omega[nu][3]
                        //                       + Wmunu[nu][3]*omega[mu][3])/2.;
                        //    // multiply term by its respective transport coefficient
                        //    term1_Vorticity = transport_coefficient4*term1_Vorticity;
                        //    // full term is
                        //    Vorticity_term = term1_Vorticity;
                        //    Vorticity_term = 0.0;
                        //} else {
                        //    Vorticity_term = 0.0;
                        //}

                        // Add nonlinear term in shear-stress tensor
                        //  transport_coefficient3*Delta(mu nu)(alpha beta)*Wmu
                        //  gamma sigma nu gamma
                        double Wsigma, Wsigma_term;
                        double term1_Wsigma, term2_Wsigma;
                        if (include_Wsigma_term == 1) {
                            Wsigma = (
                                   Wmunu[0][0]*sigma[0][0]
                                 + Wmunu[1][1]*sigma[1][1]
                                 + Wmunu[2][2]*sigma[2][2]
                                 + Wmunu[3][3]*sigma[3][3]
                                 - 2.*(  Wmunu[0][1]*sigma[0][1]
                                       + Wmunu[0][2]*sigma[0][2]
                                       + Wmunu[0][3]*sigma[0][3])
                                 +2.*(  Wmunu[1][2]*sigma[1][2]
                                      + Wmunu[1][3]*sigma[1][3]
                                      + Wmunu[2][3]*sigma[2][3]));
                            term1_Wsigma = ( - Wmunu[mu][0]*sigma[nu][0]
                                             - Wmunu[nu][0]*sigma[mu][0]
                                             + Wmunu[mu][1]*sigma[nu][1]
                                             + Wmunu[nu][1]*sigma[mu][1]
                                             + Wmunu[mu][2]*sigma[nu][2]
                                             + Wmunu[nu][2]*sigma[mu][2]
                                             + Wmunu[mu][3]*sigma[nu][3]
                                             + Wmunu[nu][3]*sigma[mu][3])/2.;

                            term2_Wsigma = (-(1./3.)*(DATA_ptr->gmunu[mu][nu]
                                                      //+ grid_pt->u[rk_flag][mu]
                                                      //  *grid_pt->u[rk_flag][nu])*Wsigma);
                                                      + vis_array[idx][15+mu]
                                                        *vis_array[idx][15+nu])
                                                     *Wsigma);
                            // multiply term by its respective transport coefficient
                            term1_Wsigma = transport_coefficient3*term1_Wsigma;
                            term2_Wsigma = transport_coefficient3*term2_Wsigma;

                            // full term is
                            Wsigma_term = -term1_Wsigma - term2_Wsigma;
                        } else {
                            Wsigma_term = 0.0;
                        }
                        // Add nonlinear term in shear-stress tensor
                        // transport_coefficient*Delta(mu nu)(alpha beta)*Wmu
                        // gamma Wnu gamma
                        double Wsquare, WW_term;
                        double term1_WW, term2_WW;
                        if (include_WWterm == 1) {
                            Wsquare = (  Wmunu[0][0]*Wmunu[0][0]
                                       + Wmunu[1][1]*Wmunu[1][1]
                                       + Wmunu[2][2]*Wmunu[2][2]
                                       + Wmunu[3][3]*Wmunu[3][3]
                                - 2.*(  Wmunu[0][1]*Wmunu[0][1]
                                      + Wmunu[0][2]*Wmunu[0][2]
                                      + Wmunu[0][3]*Wmunu[0][3])
                                + 2.*(  Wmunu[1][2]*Wmunu[1][2]
                                      + Wmunu[1][3]*Wmunu[1][3]
                                      + Wmunu[2][3]*Wmunu[2][3]));
                            term1_WW = ( - Wmunu[mu][0]*Wmunu[nu][0]
                                         + Wmunu[mu][1]*Wmunu[nu][1]
                                         + Wmunu[mu][2]*Wmunu[nu][2]
                                         + Wmunu[mu][3]*Wmunu[nu][3]);

                            term2_WW = (
                                -(1./3.)*(DATA_ptr->gmunu[mu][nu]
                                          //+ grid_pt->u[rk_flag][mu]*grid_pt->u[rk_flag][nu])
                                          + vis_array[idx][15+mu]
                                            *vis_array[idx][15+nu])
                                *Wsquare);

                            // multiply term by its respective transport coefficient
                            term1_WW = term1_WW*transport_coefficient;
                            term2_WW = term2_WW*transport_coefficient;

                            // full term is
                            // sign changes according to metric sign convention
                            WW_term = -term1_WW - term2_WW;
                        } else {
                            WW_term = 0.0;
                        }

                        // Add coupling to bulk viscous pressure
                        // transport_coefficient_b*Bulk*sigma^mu nu
                        // transport_coefficient2_b*Bulk*W^mu nu
                        double Bulk_Sigma, Bulk_Sigma_term;
                        double Bulk_W, Bulk_W_term;
                        double Coupling_to_Bulk;

                        //Bulk_Sigma = grid_pt->pi_b[rk_flag]*sigma[mu][nu];
                        //Bulk_W = grid_pt->pi_b[rk_flag]*Wmunu[mu][nu];
                        Bulk_Sigma = vis_array[idx][14]*sigma[mu][nu];
                        Bulk_W = vis_array[idx][14]*Wmunu[mu][nu];

                        // multiply term by its respective transport coefficient
                        Bulk_Sigma_term = Bulk_Sigma*transport_coefficient_b;
                        Bulk_W_term = Bulk_W*transport_coefficient2_b;

                        // full term is
                        // first term: 
                        // sign changes according to metric sign convention
                        Coupling_to_Bulk = -Bulk_Sigma_term + Bulk_W_term;

                        // final answer is
                        double SW = ((NS_term + tempf + Vorticity_term
                                      + Wsigma_term + WW_term
                                      + Coupling_to_Bulk)
                                     /(tau_pi));
                        vis_array_new[idx][idx_1d] += SW*(DATA_ptr->delta_tau);
                    }
                }
            }
        }
    }
    return(0);
}


int Diss::Make_uWRHS(double tau, int n_cell_eta, int n_cell_x, int n_cell_y,
                     double **vis_array, double **vis_nbr_x,
                     double **vis_nbr_y, double **vis_nbr_eta,
                     double **velocity_array, double **vis_array_new) {

    if (DATA_ptr->turn_on_shear == 0)
        return(1);

    int temp[6] = {4, 5, 6, 7, 8 ,9};
    vector<int> idx_1d_list(temp, temp+6);  // shear components
    if (DATA_ptr->turn_on_diff == 1) {  // diffusion components
        for (int ii = 10; ii < 14; ii++) {
            idx_1d_list.push_back(ii);
        }
    }
    if (DATA_ptr->turn_on_bulk == 1) {  // bulk viscous pressure
        idx_1d_list.push_back(14);
    }

    int mu_list[] = {1, 1, 1, 2, 2, 3};
    int nu_list[] = {1, 2, 3, 2, 3, 3};

    double Wmunu_local[4][4];
    for (int k = 0; k < n_cell_eta; k++) {
        for (int i = 0; i < n_cell_x; i++) {
            for (int j = 0; j < n_cell_y; j++) {
                int idx = j + i*n_cell_y + k*n_cell_x*n_cell_y;

                for (int aa = 0; aa < 4; aa++) {
                    for (int bb = aa; bb < 4; bb++) {
                        int idx_1d = util->map_2d_idx_to_1d(aa, bb);
                        Wmunu_local[aa][bb] = vis_array[idx][idx_1d];
                    }
                }
                for (int aa = 0; aa < 4; aa++) {
                    for (int bb = aa+1; bb < 4; bb++) {
                        Wmunu_local[bb][aa] = Wmunu_local[aa][bb];
                    }
                }

                // Kurganov-Tadmor for Wmunu */
                // implement 
                // partial_tau (utau Wmn) + (1/tau)partial_eta (ueta Wmn) 
                // + partial_x (ux Wmn) + partial_y (uy Wmn) + utau Wmn/tau
                // = SW 
                // or the right hand side of,
                // partial_tau (utau Wmn) = 
                //                  - (1/tau)partial_eta (ueta Wmn)
                //                  - partial_x (ux Wmn) - partial_y (uy Wmn) 
                //                  - utau Wmn/tau + SW*/

                // the local velocity is just u_x/u_tau, u_y/u_tau,
                //                            u_eta/tau/u_tau
                // KT flux is given by 
                // H_{j+1/2} = (fRph + fLph)/2 - ax(uRph - uLph) 
                // Here fRph = ux WmnRph and ax uRph = |ux/utau|_max utau Wmn
                // This is the second step in the operator splitting. it uses
                // rk_flag+1 as initial condition
    
                double taufactor;
                double g, gp1, gm1, gp2, gm2, a, am1, ap1, ax;
                double f, fp1, fm1, fp2, fm2;
                double uWphR, uWphL, uWmhR, uWmhL, WphR, WphL, WmhR, WmhL;
                double HWph, HWmh, HW;
                int idx_p_2, idx_p_1, idx_m_1, idx_m_2;
                for (unsigned int ii = 0; ii < idx_1d_list.size(); ii++) {
                    int idx_1d = idx_1d_list[ii];

                    // the derivative part is the same for all viscous
                    // components
                    vis_array_new[idx][idx_1d] = (
                                vis_array[idx][idx_1d]*vis_array[idx][15]);

                    double sum = 0.0;
                    // x-direction
                    taufactor = 1.0;
                    /* Get_uWmns */
                    //g = grid_pt->Wmunu[rk_flag][idx_1d];
                    //f = g*grid_pt->u[rk_flag][direc];
                    //g *= grid_pt->u[rk_flag][0];
                    g = vis_array[idx][idx_1d];
                    f = g*vis_array[idx][16];
                    g *= vis_array[idx][15];
                    //a = fabs(grid_pt->u[rk_flag][direc]);
                    a = fabs(vis_array[idx][16])/vis_array[idx][15];

                    if (i + 2 < n_cell_x) {
                        idx_p_2 = j + (i+2)*n_cell_y + k*n_cell_x*n_cell_y;
                        gp2 = vis_array[idx_p_2][idx_1d];
                        fp2 = gp2*vis_array[idx_p_2][16];
                        gp2 *= vis_array[idx_p_2][15];
                    } else {
                        idx_p_2 = 4*j + k*4*n_cell_y + 4 + i - n_cell_x;
                        gp2 = vis_nbr_x[idx_p_2][idx_1d];
                        fp2 = gp2*vis_nbr_x[idx_p_2][16];
                        gp2 *= vis_nbr_x[idx_p_2][15];
                    }

                    if (i + 1 < n_cell_x) {
                        idx_p_1 = j + (i+1)*n_cell_y + k*n_cell_x*n_cell_y;
                        gp1 = vis_array[idx_p_1][idx_1d];
                        fp1 = gp1*vis_array[idx_p_1][16];
                        gp1 *= vis_array[idx_p_1][15];
                        ap1 = (fabs(vis_array[idx_p_1][16])
                               /vis_array[idx_p_1][15]);
                    } else {
                        idx_p_1 = 4*j + k*4*n_cell_y + 2;
                        gp1 = vis_nbr_x[idx_p_1][idx_1d];
                        fp1 = gp1*vis_nbr_x[idx_p_1][16];
                        gp1 *= vis_nbr_x[idx_p_1][15];
                        ap1 = (fabs(vis_nbr_x[idx_p_1][16])
                               /vis_nbr_x[idx_p_1][15]);
                    }

                    if (i - 1 >= 0) {
                        idx_m_1 = j + (i-1)*n_cell_y + k*n_cell_x*n_cell_y;
                        gm1 = vis_array[idx_m_1][idx_1d];
                        fm1 = gm1*vis_array[idx_m_1][16];
                        gm1 *= vis_array[idx_m_1][15];
                        am1 = (fabs(vis_array[idx_m_1][16])
                               /vis_array[idx_m_1][15]);
                    } else {
                        idx_m_1 = 4*j + k*4*n_cell_y + 1;
                        gm1 = vis_nbr_x[idx_m_1][idx_1d];
                        fm1 = gm1*vis_nbr_x[idx_m_1][16];
                        gm1 *= vis_nbr_x[idx_m_1][15];
                        am1 = (fabs(vis_nbr_x[idx_m_1][16])
                               /vis_nbr_x[idx_m_1][15]);
                    }

                    if (i - 2 >= 0) {
                        idx_m_2 = j + (i-2)*n_cell_y + k*n_cell_x*n_cell_y;
                        gm2 = vis_array[idx_m_2][idx_1d];
                        fm2 = gm2*vis_array[idx_m_2][16];
                        gm2 *= vis_array[idx_m_2][15];
                    } else {
                        idx_m_2 = 4*j + k*4*n_cell_y + i;
                        gm2 = vis_nbr_x[idx_m_2][idx_1d];
                        fm2 = gm2*vis_nbr_x[idx_m_2][16];
                        gm2 *= vis_nbr_x[idx_m_2][15];
                    }
                    // MakeuWmnHalfs uWmn
                    uWphR = fp1 - 0.5*minmod->minmod_dx(fp2, fp1, f);
                    uWphL = f + 0.5*minmod->minmod_dx(fp1, f, fm1);
                    uWmhR = f - 0.5*minmod->minmod_dx(fp1, f, fm1);
                    uWmhL = fm1 + 0.5*minmod->minmod_dx(f, fm1, fm2);
                    // just Wmn
                    WphR = gp1 - 0.5*minmod->minmod_dx(gp2, gp1, g);
                    WphL = g + 0.5*minmod->minmod_dx(gp1, g, gm1);
                    WmhR = g - 0.5*minmod->minmod_dx(gp1, g, gm1);
                    WmhL = gm1 + 0.5*minmod->minmod_dx(g, gm1, gm2);

                    ax = maxi(a, ap1);
                    HWph = ((uWphR + uWphL) - ax*(WphR - WphL))*0.5;
                    ax = maxi(a, am1);
                    HWmh = ((uWmhR + uWmhL) - ax*(WmhR - WmhL))*0.5;
                    HW = (HWph - HWmh)/DATA_ptr->delta_x/taufactor;
                        
                    // make partial_i (u^i Wmn)
                    sum += -HW;
            
                    // y-direction
                    taufactor = 1.0;
                    /* Get_uWmns */
                    g = vis_array[idx][idx_1d];
                    f = g*vis_array[idx][17];
                    g *= vis_array[idx][15];
                    a = fabs(vis_array[idx][17])/vis_array[idx][15];

                    if (j + 2 < n_cell_y) {
                        idx_p_2 = j + 2 + i*n_cell_y + k*n_cell_x*n_cell_y;
                        gp2 = vis_array[idx_p_2][idx_1d];
                        fp2 = gp2*vis_array[idx_p_2][17];
                        gp2 *= vis_array[idx_p_2][15];
                    } else {
                        idx_p_2 = 4*i + 4*k*n_cell_x + 4 + j - n_cell_y;
                        gp2 = vis_nbr_y[idx_p_2][idx_1d];
                        fp2 = gp2*vis_nbr_y[idx_p_2][17];
                        gp2 *= vis_nbr_y[idx_p_2][15];
                    }

                    if (j + 1 < n_cell_y) {
                        idx_p_1 = j + 1 + i*n_cell_y + k*n_cell_x*n_cell_y;
                        gp1 = vis_array[idx_p_1][idx_1d];
                        fp1 = gp1*vis_array[idx_p_1][17];
                        gp1 *= vis_array[idx_p_1][15];
                        ap1 = (fabs(vis_array[idx_p_1][17])
                               /vis_array[idx_p_1][15]);
                    } else {
                        idx_p_1 = 4*i + 4*k*n_cell_x + 2;
                        gp1 = vis_nbr_y[idx_p_1][idx_1d];
                        fp1 = gp1*vis_nbr_y[idx_p_1][17];
                        gp1 *= vis_nbr_y[idx_p_1][15];
                        ap1 = (fabs(vis_nbr_y[idx_p_1][17])
                               /vis_nbr_y[idx_p_1][15]);
                    }
                        
                    if (j - 1 >= 0) {
                        idx_m_1 = j - 1 + i*n_cell_y + k*n_cell_x*n_cell_y;
                        gm1 = vis_array[idx_m_1][idx_1d];
                        fm1 = gm1*vis_array[idx_m_1][17];
                        gm1 *= vis_array[idx_m_1][15];
                        am1 = (fabs(vis_array[idx_m_1][17])
                               /vis_array[idx_m_1][15]);
                    } else {
                        idx_m_1 = 4*i + 4*k*n_cell_x + 1;
                        gm1 = vis_nbr_y[idx_m_1][idx_1d];
                        fm1 = gm1*vis_nbr_y[idx_m_1][17];
                        gm1 *= vis_nbr_y[idx_m_1][15];
                        am1 = (fabs(vis_nbr_y[idx_m_1][17])
                               /vis_nbr_y[idx_m_1][15]);
                    }

                    if (j - 2 >= 0) {
                        idx_m_2 = j - 2 + i*n_cell_y + k*n_cell_x*n_cell_y;
                        gm2 = vis_array[idx_m_2][idx_1d];
                        fm2 = gm2*vis_array[idx_m_2][17];
                        gm2 *= vis_array[idx_m_2][15];
                    } else {
                        idx_m_2 = 4*i + 4*k*n_cell_x + j;
                        gm2 = vis_nbr_y[idx_m_2][idx_1d];
                        fm2 = gm2*vis_nbr_y[idx_m_2][17];
                        gm2 *= vis_nbr_y[idx_m_2][15];
                    }
                    // MakeuWmnHalfs uWmn
                    uWphR = fp1 - 0.5*minmod->minmod_dx(fp2, fp1, f);
                    uWphL = f + 0.5*minmod->minmod_dx(fp1, f, fm1);
                    uWmhR = f - 0.5*minmod->minmod_dx(fp1, f, fm1);
                    uWmhL = fm1 + 0.5*minmod->minmod_dx(f, fm1, fm2);
                    // just Wmn
                    WphR = gp1 - 0.5*minmod->minmod_dx(gp2, gp1, g);
                    WphL = g + 0.5*minmod->minmod_dx(gp1, g, gm1);
                    WmhR = g - 0.5*minmod->minmod_dx(gp1, g, gm1);
                    WmhL = gm1 + 0.5*minmod->minmod_dx(g, gm1, gm2);
                    ax = maxi(a, ap1);
                    HWph = ((uWphR + uWphL) - ax*(WphR - WphL))*0.5;
                    ax = maxi(a, am1);
                    HWmh = ((uWmhR + uWmhL) - ax*(WmhR - WmhL))*0.5;
                    HW = (HWph - HWmh)/DATA_ptr->delta_y/taufactor;
                    // make partial_i (u^i Wmn)
                    sum += -HW;
            
                    // eta-direction
                    taufactor = tau;
                    /* Get_uWmns */
                    g = vis_array[idx][idx_1d];
                    f = g*vis_array[idx][18];
                    g *= vis_array[idx][15];
                    a = fabs(vis_array[idx][18])/vis_array[idx][15];

                    if (k + 2 < n_cell_eta) {
                        idx_p_2 = j + i*n_cell_y + (k+2)*n_cell_x*n_cell_y;
                        gp2 = vis_array[idx_p_2][idx_1d];
                        fp2 = gp2*vis_array[idx_p_2][18];
                        gp2 *= vis_array[idx_p_2][15];
                    } else {
                        idx_p_2 = 4*j + 4*i*n_cell_y + 4 + k - n_cell_eta;
                        gp2 = vis_nbr_eta[idx_p_2][idx_1d];
                        fp2 = gp2*vis_nbr_eta[idx_p_2][18];
                        gp2 *= vis_nbr_eta[idx_p_2][15];
                    }

                    if (k + 1 < n_cell_eta) {
                        idx_p_1 = j + i*n_cell_y + (k+1)*n_cell_x*n_cell_y;
                        gp1 = vis_array[idx_p_1][idx_1d];
                        fp1 = gp1*vis_array[idx_p_1][18];
                        gp1 *= vis_array[idx_p_1][15];
                        ap1 = (fabs(vis_array[idx_p_1][18])
                                /vis_array[idx_p_1][15]);
                    } else {
                        idx_p_1 = 4*j + 4*i*n_cell_y + 2;
                        gp1 = vis_nbr_eta[idx_p_1][idx_1d];
                        fp1 = gp1*vis_nbr_eta[idx_p_1][18];
                        gp1 *= vis_nbr_eta[idx_p_1][15];
                        ap1 = (fabs(vis_nbr_eta[idx_p_1][18])
                                /vis_nbr_eta[idx_p_1][15]);
                    }

                    if (k - 1 >= 0) {
                        idx_m_1 = j + i*n_cell_y + (k-1)*n_cell_x*n_cell_y;
                        gm1 = vis_array[idx_m_1][idx_1d];
                        fm1 = gm1*vis_array[idx_m_1][18];
                        gm1 *= vis_array[idx_m_1][15];
                        am1 = (fabs(vis_array[idx_m_1][18])
                                /vis_array[idx_m_1][15]);
                    } else {
                        idx_m_1 = 4*j + 4*i*n_cell_y + 1;
                        gm1 = vis_nbr_eta[idx_m_1][idx_1d];
                        fm1 = gm1*vis_nbr_eta[idx_m_1][18];
                        gm1 *= vis_nbr_eta[idx_m_1][15];
                        am1 = (fabs(vis_nbr_eta[idx_m_1][18])
                                /vis_nbr_eta[idx_m_1][15]);
                    }

                    if (k - 2 >= 0) {
                        idx_m_2 = j + i*n_cell_y + (k-2)*n_cell_x*n_cell_y;
                        gm2 = vis_array[idx_m_2][idx_1d];
                        fm2 = gm2*vis_array[idx_m_2][18];
                        gm2 *= vis_array[idx_m_2][15];
                    } else {
                        idx_m_2 = 4*j + 4*i*n_cell_y + k;
                        gm2 = vis_nbr_eta[idx_m_2][idx_1d];
                        fm2 = gm2*vis_nbr_eta[idx_m_2][18];
                        gm2 *= vis_nbr_eta[idx_m_2][15];
                    }
                    // MakeuWmnHalfs uWmn
                    uWphR = fp1 - 0.5*minmod->minmod_dx(fp2, fp1, f);
                    uWphL = f + 0.5*minmod->minmod_dx(fp1, f, fm1);
                    uWmhR = f - 0.5*minmod->minmod_dx(fp1, f, fm1);
                    uWmhL = fm1 + 0.5*minmod->minmod_dx(f, fm1, fm2);
                    // just Wmn
                    WphR = gp1 - 0.5*minmod->minmod_dx(gp2, gp1, g);
                    WphL = g + 0.5*minmod->minmod_dx(gp1, g, gm1);
                    WmhR = g - 0.5*minmod->minmod_dx(gp1, g, gm1);
                    WmhL = gm1 + 0.5*minmod->minmod_dx(g, gm1, gm2);
                    ax = maxi(a, ap1);
                    HWph = ((uWphR + uWphL) - ax*(WphR - WphL))*0.5;
                    ax = maxi(a, am1);
                    HWmh = ((uWmhR + uWmhL) - ax*(WmhR - WmhL))*0.5;
                    HW = (HWph - HWmh)/DATA_ptr->delta_eta/taufactor;
                    // make partial_i (u^i Wmn)
                    sum += -HW;

                    // the following geometric parts are different for
                    // individual pi^\mu\nu, q^\mu, and Pi

                    if (idx_1d < 10) {
                        // geometric terms for shear pi^\mu\nu
                        int mu = mu_list[ii];
                        int nu = nu_list[ii];
                        // add a source term -u^tau Wmn/tau
                        //   due to the coordinate change to tau-eta
                        //sum += (- (grid_pt->u[rk_flag][0]*Wmunu_local[mu][nu])/tau
                        //        + (theta_local*Wmunu_local[mu][nu]));
                        sum += (- (vis_array[idx][15]*Wmunu_local[mu][nu])/tau
                                + (velocity_array[idx][0]*Wmunu_local[mu][nu])
                               );

                        // this is from udW = d(uW) - Wdu = RHS
                        // or d(uW) = udW + Wdu
                        // this term is being added to the rhs so that
                        // -4/3 + 1 = -1/3
                        // other source terms due to the coordinate
                        // change to tau-eta
                        double tempf = 0.0;
                        tempf = (
                            - (DATA_ptr->gmunu[3][mu])*(Wmunu_local[0][nu])
                            - (DATA_ptr->gmunu[3][nu])*(Wmunu_local[0][mu])
                            + (DATA_ptr->gmunu[0][mu])*(Wmunu_local[3][nu])
                            + (DATA_ptr->gmunu[0][nu])*(Wmunu_local[3][mu])
                            + (Wmunu_local[3][nu])
                              //*(grid_pt->u[rk_flag][mu])*(grid_pt->u[rk_flag][0])
                              *(vis_array[idx][15+mu])*(vis_array[idx][15])
                            + (Wmunu_local[3][mu])
                              //*(grid_pt->u[rk_flag][nu])*(grid_pt->u[rk_flag][0])
                              *(vis_array[idx][15+nu])*(vis_array[idx][15])
                            - (Wmunu_local[0][nu])
                              //*(grid_pt->u[rk_flag][mu])*(grid_pt->u[rk_flag][3])
                              *(vis_array[idx][15+mu])*(vis_array[idx][18])
                            - (Wmunu_local[0][mu])
                              //*(grid_pt->u[rk_flag][nu])*(grid_pt->u[rk_flag][3]))
                              //*(grid_pt->u[rk_flag][3]/tau);
                              *(vis_array[idx][15+nu])*(vis_array[idx][18]))
                              *(vis_array[idx][18]/tau);
                        for (int ic = 0; ic < 4; ic++) {
                            double ic_fac = (ic == 0 ? -1.0 : 1.0);
                            //tempf += (
                            //      (Wmunu_local[ic][nu])*(grid_pt->u[rk_flag][mu])
                            //       *(a_local[ic])*ic_fac
                            //    + (Wmunu_local[ic][mu])*(grid_pt->u[rk_flag][nu])
                            //       *(a_local[ic])*ic_fac);
                            tempf += (
                                  (Wmunu_local[ic][nu])*(vis_array[idx][15+mu])
                                   *(velocity_array[idx][1+ic])*ic_fac
                                + (Wmunu_local[ic][mu])*(vis_array[idx][15+nu])
                                   *(velocity_array[idx][1+ic])*ic_fac);
                        }
                        sum += tempf;
                    } else if (idx_1d == 14) {
                        // geometric terms for bulk Pi
                        //sum -= (pi_b[rk_flag])*(u[rk_flag][0])/tau;
                        //sum += (pi_b[rk_flag])*theta_local;
                        sum -= vis_array[idx][14]*vis_array[idx][15]/tau;
                        sum += vis_array[idx][14]*velocity_array[idx][0];
                    }
                    //w_rhs[mu][nu] = sum*(DATA_ptr->delta_tau);
                    vis_array_new[idx][idx_1d] += (
                                        sum*(DATA_ptr->delta_tau));
                }
            }
        }
    }
    return(1);
}


int Diss::Make_uPRHS(double tau, double *p_rhs,
                     double **vis_array, double **vis_nbr_x,
                     double **vis_nbr_y, double **vis_nbr_eta,
                     double **velocity_array) {

    /* Kurganov-Tadmor for Pi */
    /* implement 
      partial_tau (utau Pi) + (1/tau)partial_eta (ueta Pi) 
      + partial_x (ux Pi) + partial_y (uy Pi) + utau Pi/tau = SP 
      or the right hand side of
      partial_tau (utau Pi) = -
      (1/tau)partial_eta (ueta Pi) - partial_x (ux Pi) - partial_y (uy Pi)
      - utau Pi/tau + SP 
      */
    
    /* the local velocity is just u_x/u_tau, u_y/u_tau, u_eta/tau/u_tau */
    /* KT flux is given by 
       H_{j+1/2} = (fRph + fLph)/2 - ax(uRph - uLph) 
       Here fRph = ux PiRph and ax uRph = |ux/utau|_max utau Pin */
    
    /* This is the second step in the operator splitting. it uses
       rk_flag+1 as initial condition */

    if (DATA_ptr->turn_on_bulk == 0) {
        *p_rhs = 0.0;
        return(0);
    }

    int idx = 0;

    double f, fp1, fm1, fp2, fm2;
    double g, gp1, gm1, gp2, gm2, a, am1, ap1, ax;
    double uPiphR, uPiphL, uPimhR, uPimhL, PiphR, PiphL, PimhR, PimhL;
    double HPiph, HPimh, taufactor, HPi;
    double sum = 0.0;

    // x-direction
    taufactor = 1.0;
    // Get_uPis
    //g = grid_pt->pi_b[rk_flag];
    //f = g*grid_pt->u[rk_flag][direc];
    //g *= grid_pt->u[rk_flag][0];
    g = vis_array[idx][14];
    f = g*vis_array[idx][16];
    g *= vis_array[idx][15];

    //gp2 = grid_pt->nbr_p_2[direc]->pi_b[rk_flag];
    //fp2 = gp2*grid_pt->nbr_p_2[direc]->u[rk_flag][direc];
    //gp2 *= grid_pt->nbr_p_2[direc]->u[rk_flag][0];
    gp2 = vis_nbr_x[3][14];
    fp2 = gp2*vis_nbr_x[3][16];
    gp2 *= vis_nbr_x[3][15];
    
    //gp1 = grid_pt->nbr_p_1[direc]->pi_b[rk_flag];
    //fp1 = gp1*grid_pt->nbr_p_1[direc]->u[rk_flag][direc];
    //gp1 *= grid_pt->nbr_p_1[direc]->u[rk_flag][0];
    gp1 = vis_nbr_x[2][14];
    fp1 = gp1*vis_nbr_x[2][16];
    gp1 *= vis_nbr_x[2][15];
    
    //gm1 = grid_pt->nbr_m_1[direc]->pi_b[rk_flag];
    //fm1 = gm1*grid_pt->nbr_m_1[direc]->u[rk_flag][direc];
    //gm1 *= grid_pt->nbr_m_1[direc]->u[rk_flag][0];
    gm1 = vis_nbr_x[1][14];
    fm1 = gm1*vis_nbr_x[1][16];
    gm1 *= vis_nbr_x[1][15];
    
    //gm2 = grid_pt->nbr_m_2[direc]->pi_b[rk_flag];
    //fm2 = gm2*grid_pt->nbr_m_2[direc]->u[rk_flag][direc];
    //gm2 *= grid_pt->nbr_m_2[direc]->u[rk_flag][0];
    gm2 = vis_nbr_x[0][14];
    fm2 = gm2*vis_nbr_x[0][16];
    gm2 *= vis_nbr_x[0][15];

    //  Make upi Halfs uPi
    uPiphR = fp1 - 0.5*minmod->minmod_dx(fp2, fp1, f); 
    uPiphL = f + 0.5*minmod->minmod_dx(fp1, f, fm1);
    uPimhR = f - 0.5*minmod->minmod_dx(fp1, f, fm1);
    uPimhL = fm1 + 0.5*minmod->minmod_dx(f, fm1, fm2);

    // just Pi
    PiphR = gp1 - 0.5*minmod->minmod_dx(gp2, gp1, g); 
    PiphL = g + 0.5*minmod->minmod_dx(gp1, g, gm1);
    PimhR = g - 0.5*minmod->minmod_dx(gp1, g, gm1);
    PimhL = gm1 + 0.5*minmod->minmod_dx(g, gm1, gm2);

    // MakePimnCurrents following Kurganov-Tadmor
    //a = fabs(grid_pt->u[rk_flag][direc]);
    //a /= grid_pt->u[rk_flag][0];
    a = fabs(vis_array[idx][16])/vis_array[idx][15];
    //am1 = fabs(grid_pt->nbr_m_1[direc]->u[rk_flag][direc]);
    //am1 /= grid_pt->nbr_m_1[direc]->u[rk_flag][0];
    am1 = fabs(vis_nbr_x[1][16])/vis_nbr_x[1][15];
    //ap1 = fabs(grid_pt->nbr_p_1[direc]->u[rk_flag][direc]);
    //ap1 /= grid_pt->nbr_p_1[direc]->u[rk_flag][0];
    ap1 = fabs(vis_nbr_x[2][16])/vis_nbr_x[2][15];
    
    ax = maxi(a, ap1);
    HPiph = ((uPiphR + uPiphL) - ax*(PiphR - PiphL))*0.5;
    ax = maxi(a, am1); 
    HPimh = ((uPimhR + uPimhL) - ax*(PimhR - PimhL))*0.5;
    HPi = (HPiph - HPimh)/DATA_ptr->delta_x/taufactor;
    // make partial_i (u^i Pi)
    sum += -HPi;

    // y-direction
    taufactor = 1.0;
    // Get_uPis
    g = vis_array[idx][14];
    f = g*vis_array[idx][17];
    g *= vis_array[idx][15];

    gp2 = vis_nbr_y[3][14];
    fp2 = gp2*vis_nbr_y[3][17];
    gp2 *= vis_nbr_y[3][15];
    
    gp1 = vis_nbr_y[2][14];
    fp1 = gp1*vis_nbr_y[2][17];
    gp1 *= vis_nbr_y[2][15];
    
    gm1 = vis_nbr_y[1][14];
    fm1 = gm1*vis_nbr_y[1][17];
    gm1 *= vis_nbr_y[1][15];
    
    gm2 = vis_nbr_y[0][14];
    fm2 = gm2*vis_nbr_y[0][17];
    gm2 *= vis_nbr_y[0][15];

    //  Make upi Halfs uPi
    uPiphR = fp1 - 0.5*minmod->minmod_dx(fp2, fp1, f); 
    uPiphL = f + 0.5*minmod->minmod_dx(fp1, f, fm1);
    uPimhR = f - 0.5*minmod->minmod_dx(fp1, f, fm1);
    uPimhL = fm1 + 0.5*minmod->minmod_dx(f, fm1, fm2);

    // just Pi
    PiphR = gp1 - 0.5*minmod->minmod_dx(gp2, gp1, g); 
    PiphL = g + 0.5*minmod->minmod_dx(gp1, g, gm1);
    PimhR = g - 0.5*minmod->minmod_dx(gp1, g, gm1);
    PimhL = gm1 + 0.5*minmod->minmod_dx(g, gm1, gm2);

    // MakePimnCurrents following Kurganov-Tadmor
    a = fabs(vis_array[idx][17])/vis_array[idx][15];
    am1 = fabs(vis_nbr_y[1][17])/vis_nbr_y[1][15];
    ap1 = fabs(vis_nbr_y[2][17])/vis_nbr_y[2][15];
    ax = maxi(a, ap1);
    HPiph = ((uPiphR + uPiphL) - ax*(PiphR - PiphL))*0.5;
    ax = maxi(a, am1); 
    HPimh = ((uPimhR + uPimhL) - ax*(PimhR - PimhL))*0.5;
    HPi = (HPiph - HPimh)/DATA_ptr->delta_y/taufactor;
    // make partial_i (u^i Pi)
    sum += -HPi;
    
    // eta-direction
    taufactor = tau;
    // Get_uPis
    g = vis_array[idx][14];
    f = g*vis_array[idx][18];
    g *= vis_array[idx][15];

    gp2 = vis_nbr_eta[3][14];
    fp2 = gp2*vis_nbr_eta[3][18];
    gp2 *= vis_nbr_eta[3][15];
    
    gp1 = vis_nbr_eta[2][14];
    fp1 = gp1*vis_nbr_eta[2][18];
    gp1 *= vis_nbr_eta[2][15];
    
    gm1 = vis_nbr_eta[1][14];
    fm1 = gm1*vis_nbr_eta[1][18];
    gm1 *= vis_nbr_eta[1][15];
    
    gm2 = vis_nbr_eta[0][14];
    fm2 = gm2*vis_nbr_eta[0][18];
    gm2 *= vis_nbr_eta[0][15];

    //  Make upi Halfs uPi
    uPiphR = fp1 - 0.5*minmod->minmod_dx(fp2, fp1, f); 
    uPiphL = f + 0.5*minmod->minmod_dx(fp1, f, fm1);
    uPimhR = f - 0.5*minmod->minmod_dx(fp1, f, fm1);
    uPimhL = fm1 + 0.5*minmod->minmod_dx(f, fm1, fm2);

    // just Pi
    PiphR = gp1 - 0.5*minmod->minmod_dx(gp2, gp1, g); 
    PiphL = g + 0.5*minmod->minmod_dx(gp1, g, gm1);
    PimhR = g - 0.5*minmod->minmod_dx(gp1, g, gm1);
    PimhL = gm1 + 0.5*minmod->minmod_dx(g, gm1, gm2);

    // MakePimnCurrents following Kurganov-Tadmor
    a = fabs(vis_array[idx][18])/vis_array[idx][15];
    am1 = fabs(vis_nbr_eta[1][18])/vis_nbr_eta[1][15];
    ap1 = fabs(vis_nbr_eta[2][18])/vis_nbr_eta[2][15];
    ax = maxi(a, ap1);
    HPiph = ((uPiphR + uPiphL) - ax*(PiphR - PiphL))*0.5;
    ax = maxi(a, am1); 
    HPimh = ((uPimhR + uPimhL) - ax*(PimhR - PimhL))*0.5;
    HPi = (HPiph - HPimh)/DATA_ptr->delta_eta/taufactor;
    // make partial_i (u^i Pi)
    sum += -HPi;
    
    // add a source term due to the coordinate change to tau-eta
    //sum -= (grid_pt->pi_b[rk_flag])*(grid_pt->u[rk_flag][0])/tau;
    //sum += (grid_pt->pi_b[rk_flag])*theta_local;
    sum -= vis_array[idx][14]*vis_array[idx][15]/tau;
    sum += vis_array[idx][14]*velocity_array[idx][0];
    *p_rhs = sum*(DATA_ptr->delta_tau);
    return(0);
}


double Diss::Make_uPiSource(double tau, int n_cell_eta, int n_cell_x,
                            int n_cell_y,
                            double **vis_array, double **velocity_array,
                            double **grid_array, double **vis_array_new) {
    // switch to include non-linear coupling terms in the bulk pi evolution
    int include_BBterm = 1;
    int include_coupling_to_shear = 1;
 
    if (DATA_ptr->turn_on_bulk == 0) return 0.0;

    for (int k = 0; k < n_cell_eta; k++) {
        for (int i = 0; i < n_cell_x; i++) {
            for (int j = 0; j < n_cell_y; j++) {
                int idx = j + i*n_cell_y + k*n_cell_x*n_cell_y;
                // defining bulk viscosity coefficient
                double epsilon = grid_array[idx][0];
                double rhob = grid_array[idx][4];
                double temperature = eos->get_temperature(epsilon, rhob);

                // cs2 is the velocity of sound squared
                double cs2 = eos->get_cs2(epsilon, rhob);  
                double pressure = eos->get_pressure(epsilon, rhob);

                // T dependent bulk viscosity from Gabriel
                double bulk =
                            get_temperature_dependent_zeta_s(temperature);
                bulk = bulk*(epsilon + pressure)/temperature;

                // defining bulk relaxation time and
                // additional transport coefficients
                // Bulk relaxation time from kinetic theory
                double Bulk_Relax_time = (1./14.55/(1./3.-cs2)/(1./3.-cs2)
                                          /(epsilon + pressure)*bulk);

                // from kinetic theory, small mass limit
                double transport_coeff1 = 2.0/3.0*(Bulk_Relax_time);
                double transport_coeff2 = 0.;  // not known; put 0

                // from kinetic theory
                double transport_coeff1_s = (8./5.*(1./3.-cs2)
                                             *Bulk_Relax_time);
                double transport_coeff2_s = 0.;  // not known;  put 0

                // Computing Navier-Stokes term (-bulk viscosity * theta)
                //double NS_term = -bulk*theta_local;
                double NS_term = -bulk*velocity_array[idx][0];

                // Computing relaxation term and nonlinear term:
                // - Bulk - transport_coeff1*Bulk*theta
                //double tempf = (-(grid_pt->pi_b[rk_flag])
                //         - transport_coeff1*theta_local
                //           *(grid_pt->pi_b[rk_flag]));
                double tempf = (- vis_array[idx][14]
                    - transport_coeff1*velocity_array[idx][0]
                      *vis_array[idx][14]);

                // Computing nonlinear term: + transport_coeff2*Bulk*Bulk
                double BB_term = 0.0;
                if (include_BBterm == 1) {
                    //BB_term = (transport_coeff2*(grid_pt->pi_b[rk_flag])
                    //           *(grid_pt->pi_b[rk_flag]));
                    BB_term = (transport_coeff2*vis_array[idx][14]
                               *vis_array[idx][14]);
                }

                // Computing terms that Couple with shear-stress tensor
                double Wsigma, WW, Shear_Sigma_term, Shear_Shear_term;
                double Coupling_to_Shear;

                if (include_coupling_to_shear == 1) {
                    // Computing sigma^mu^nu
                    double sigma[4][4], Wmunu[4][4];
                    for (int a = 0; a < 4 ; a++) {
                        for (int b = a; b < 4; b++) {
                            int idx_1d = util->map_2d_idx_to_1d(a, b);
                            sigma[a][b] = velocity_array[idx][6+idx_1d];
                            Wmunu[a][b] = vis_array[idx][idx_1d];
                        }
                    }

                    Wsigma = (  Wmunu[0][0]*sigma[0][0]
                              + Wmunu[1][1]*sigma[1][1]
                              + Wmunu[2][2]*sigma[2][2]
                              + Wmunu[3][3]*sigma[3][3]
                              - 2.*(  Wmunu[0][1]*sigma[0][1]
                                    + Wmunu[0][2]*sigma[0][2]
                                    + Wmunu[0][3]*sigma[0][3])
                              + 2.*(  Wmunu[1][2]*sigma[1][2]
                                    + Wmunu[1][3]*sigma[1][3]
                                    + Wmunu[2][3]*sigma[2][3]));

                    WW = (   Wmunu[0][0]*Wmunu[0][0]
                           + Wmunu[1][1]*Wmunu[1][1]
                           + Wmunu[2][2]*Wmunu[2][2]
                           + Wmunu[3][3]*Wmunu[3][3]
                           - 2.*(  Wmunu[0][1]*Wmunu[0][1]
                                 + Wmunu[0][2]*Wmunu[0][2]
                                 + Wmunu[0][3]*Wmunu[0][3])
                           + 2.*(  Wmunu[1][2]*Wmunu[1][2]
                                 + Wmunu[1][3]*Wmunu[1][3]
                                 + Wmunu[2][3]*Wmunu[2][3]));
                    // multiply term by its transport coefficient
                    Shear_Sigma_term = Wsigma*transport_coeff1_s;
                    Shear_Shear_term = WW*transport_coeff2_s;

                    // full term that couples to shear is
                    Coupling_to_Shear = (- Shear_Sigma_term
                                         + Shear_Shear_term);
                } else {
                    Coupling_to_Shear = 0.0;
                }
                
                // Final Answer
                double Final_Answer = (NS_term + tempf + BB_term
                                        + Coupling_to_Shear)/Bulk_Relax_time;
                vis_array_new[idx][14] += Final_Answer*(DATA_ptr->delta_tau);
            }
        }
    }
    return(0);
}


/* Sangyong Nov 18 2014 */
/* baryon current parts */
/* this contains the source terms
   that is, all the terms that are not part of the current */
/* for the q part, we don't do tau*u*q we just do u*q 
   this part contains 
    -(1/tau_rho)(q[a] + kappa g[a][b]Dtildemu[b]
                 + kappa u[a] u[b]g[b][c]Dtildemu[c])
    +Delta[a][tau] u[eta] q[eta]/tau
    -Delta[a][eta] u[eta] q[tau]/tau
    -u[a]u[b]g[b][e] Dq[e]
*/
double Diss::Make_uqSource(double tau, int n_cell_eta, int n_cell_x,
                           int n_cell_y,
                           double **vis_array, double **velocity_array,
                           double **grid_array, double **vis_array_new) {

    if (DATA_ptr->turn_on_diff == 0) return 0.0;
 
    double sigma[4][4];

    for (int k = 0; k < n_cell_eta; k++) {
        for (int i = 0; i < n_cell_x; i++) {
            for (int j = 0; j < n_cell_y; j++) {
                int idx = j + i*n_cell_y + k*n_cell_x*n_cell_y;

                // Useful variables to define
                double epsilon = grid_array[idx][0];
                double rhob = grid_array[idx][4];
                double pressure = eos->get_pressure(epsilon, rhob);
                double T = eos->get_temperature(epsilon, rhob);

                double kappa_coefficient = DATA_ptr->kappa_coefficient;
                double tau_rho = kappa_coefficient/(T + 1e-15);
                double mub = eos->get_mu(epsilon, rhob);
                double alpha = mub/T;
                double kappa = (kappa_coefficient
                                *(rhob/(3.*T*tanh(alpha) + 1e-15)
                                  - rhob*rhob/(epsilon + pressure)));


                for (int ii = 0; ii < 4; ii++) {
                    for (int jj = ii; jj < 4; jj++) {
                        int idx_1d = util->map_2d_idx_to_1d(ii, jj);
                        //sigma[ii][jj] = sigma_1d[idx_1d];
                        sigma[ii][jj] = velocity_array[idx][6+idx_1d];
                    }
                }
                for (int ii = 0; ii < 4; ii++) {
                    for (int jj = ii+1; jj < 4; jj++) {
                        sigma[jj][ii] = sigma[ii][jj];
                    }
                }

                // add a new non-linear term (- q \theta)
                // from conformal kinetic theory
                double transport_coeff = 1.0*tau_rho;
                // add a new non-linear term (-q^\mu \sigma_\mu\nu)
                // from 14-momentum massless
                double transport_coeff_2 = 3./5.*tau_rho;

                /* -(1/tau_rho)(q[a] + kappa g[a][b]Dtildemu[b] 
                 *              + kappa u[a] u[b]g[b][c]Dtildemu[c])
                 * + theta q[a] - q[a] u^\tau/tau
                 * + Delta[a][tau] u[eta] q[eta]/tau
                 * - Delta[a][eta] u[eta] q[tau]/tau
                 * - u[a] u[b]g[b][e] Dq[e] -> u[a] q[e] g[e][b] Du[b]
                */    
 
                // first: (1/tau_rho) part
                // recall that dUsup[4][i] = partial_i (muB/T) 
                // and dUsup[4][0] = -partial_tau (muB/T) = partial^tau (muB/T)
                // and a[4] = u^a partial_a (muB/T) = DmuB/T
                // -(1/tau_rho)(q[a] + kappa g[a][b]DmuB/T[b] 
                // + kappa u[a] u[b]g[b][c]DmuB/T[c])
                // a = nu 
                //double NS = kappa*(grid_pt->dUsup[0][4][nu] 
                //                       + grid_pt->u[rk_flag][nu]*a_local[4]);
                for (int nu = 1; nu < 4; nu++) {
                    int idx_1d = util->map_2d_idx_to_1d(4, nu);
                    //double NS = kappa*(dUsup[0][4][nu] 
                    //                       + u[rk_flag][nu]*a_local[4]);
                    double NS = kappa*(velocity_array[idx][16+nu]
                           + vis_array[idx][15+nu]*velocity_array[idx][5]);
                    if (isnan(NS)) {
                        cout << "Navier Stock term is nan! " << endl;
                        cout << vis_array[idx][idx_1d] << endl;
                        // derivative already upper index
                        cout << velocity_array[idx][16+nu] << endl;
                        cout << velocity_array[5] << endl;
                        cout << tau_rho << endl;
                        cout << kappa << endl;
                        cout << vis_array[idx][15+nu] << endl;
                    }
  
                    //double Nonlinear1 = -transport_coeff*q[nu]*theta_local;
                    double Nonlinear1 = (- transport_coeff
                                           *vis_array[idx][idx_1d]
                                           *velocity_array[idx][0]);

                    double temptemp = 0.0;
                    for (int iii = 0 ; iii < 4; iii++) {
                        temptemp += (vis_array[idx][10+iii]*sigma[iii][nu]
                                     *DATA_ptr->gmunu[iii][iii]);
                    }
                    double Nonlinear2 = - transport_coeff_2*temptemp;

                    double SW = ((-vis_array[idx][idx_1d]
                                  - NS + Nonlinear1 + Nonlinear2)
                                 /(tau_rho + 1e-15));

                    // for 1+1D numerical test
                    // SW = (-q[nu] - NS)/(tau_rho + 1e-15);

                    // all other geometric terms....
                    // + theta q[a] - q[a] u^\tau/tau
                    //SW += (theta_local - grid_pt->u[rk_flag][0]/tau)*q[nu];
                    SW += (velocity_array[idx][0]
                            - vis_array[idx][15]/tau)*vis_array[idx][idx_1d];
 
                    if (isnan(SW)) {
                        cout << "theta term is nan! " << endl;
                    }

                    // +Delta[a][tau] u[eta] q[eta]/tau 
                    double tempf = ((DATA_ptr->gmunu[nu][0] 
                                    //+ grid_pt->u[rk_flag][nu]*grid_pt->u[rk_flag][0])
                                    //  *grid_pt->u[rk_flag][3]*q[3]/tau
                                    + vis_array[idx][15+nu]*vis_array[idx][15])
                                      *vis_array[idx][18]
                                      *vis_array[idx][13]/tau
                                    - (DATA_ptr->gmunu[nu][3]
                                      // + grid_pt->u[rk_flag][nu]*grid_pt->u[rk_flag][3])
                                      //*grid_pt->u[rk_flag][3]*q[0]/tau);
                                       + vis_array[idx][15+nu]
                                         *vis_array[idx][18])
                                      *vis_array[idx][18]
                                      *vis_array[idx][10]/tau);
                    SW += tempf;
 
                    if (isnan(tempf)) {
                        cout << "Delta^{a \tau} and Delta^{a \eta} terms "
                             << "are nan!" << endl;
                    }

                    // -u[a] u[b]g[b][e] Dq[e] -> u[a] (q[e] g[e][b] Du[b])
                    tempf = 0.0;
                    for (int iii = 0; iii < 4; iii++) {
                        //tempf += q[i]*gmn(i)*a_local[i];
                        tempf += (vis_array[idx][10+iii]*gmn(iii)
                                  *velocity_array[idx][1+iii]);
                    }
                    //SW += (grid_pt->u[rk_flag][nu])*tempf;
                    SW += vis_array[idx][15+nu]*tempf;
                    
                    if (isnan(tempf)) {
                        cout << "u^a q_b Du^b term is nan! " << endl;
                    }
                    vis_array_new[idx][idx_1d] += SW*(DATA_ptr->delta_tau);
                }
            }
        }
    }
    return(0);
}


int Diss::Make_uqRHS(double tau, double **w_rhs, double **vis_array,
                     double **vis_nbr_x, double **vis_nbr_y,
                     double **vis_nbr_eta) {

    /* Kurganov-Tadmor for q */
    /* implement 
      partial_tau (utau qmu) + (1/tau)partial_eta (ueta qmu) 
      + partial_x (ux qmu) + partial_y (uy qmu) + utau qmu/tau = SW 
    or the right hand side of,
      partial_tau (utau qmu) = 
      - (1/tau)partial_eta (ueta qmu) - partial_x (ux qmu) - partial_y (uy qmu) 
      - utau qmu/tau 
    */

    /* the local velocity is just u_x/u_tau, u_y/u_tau, u_eta/tau/u_tau */
    /* KT flux is given by 
       H_{j+1/2} = (fRph + fLph)/2 - ax(uRph - uLph) 
       Here fRph = ux WmnRph and ax uRph = |ux/utau|_max utau Wmn */

    /* This is the second step in the operator splitting. it uses
       rk_flag+1 as initial condition */

    double f, fp1, fm1, fp2, fm2;
    double g, gp1, gm1, gp2, gm2, a, am1, ap1, ax;
    double uWphR, uWphL, uWmhR, uWmhL, WphR, WphL, WmhR, WmhL;
    double HWph, HWmh, HW;

    int idx = 0;
    // we use the Wmunu[4][nu] = q[nu] 
    int mu = 4;
    double taufactor = 1.0;
    for (int nu = 1; nu < 4; nu++) {
        int idx_1d = util->map_2d_idx_to_1d(mu, nu);
        double sum = 0.0;

        // x-direction
        taufactor = 1.0;
        /* Get_uWmns */
        g = vis_array[idx][idx_1d];
        f = g*vis_array[idx][16];
        g *= vis_array[idx][15];
           
        gp2 = vis_nbr_x[3][idx_1d];
        fp2 = gp2*vis_nbr_x[3][16];
        gp2 *= vis_nbr_x[3][15];
        
        gp1 = vis_nbr_x[2][idx_1d];
        fp1 = gp1*vis_nbr_x[2][16];
        gp1 *= vis_nbr_x[2][15];
        
        gm1 = vis_nbr_x[1][idx_1d];
        fm1 = gm1*vis_nbr_x[1][16];
        gm1 *= vis_nbr_x[1][15];
        
        gm2 = vis_nbr_x[0][idx_1d];
        fm2 = gm2*vis_nbr_x[0][16];
        gm2 *= vis_nbr_x[0][15];

        // MakeuWmnHalfs uWmn
        uWphR = fp1 - 0.5*minmod->minmod_dx(fp2, fp1, f);
        uWphL = f + 0.5*minmod->minmod_dx(fp1, f, fm1);
        uWmhR = f - 0.5*minmod->minmod_dx(fp1, f, fm1);
        uWmhL = fm1 + 0.5*minmod->minmod_dx(f, fm1, fm2);

        // just Wmn
        WphR = gp1 - 0.5*minmod->minmod_dx(gp2, gp1, g);
        WphL = g + 0.5*minmod->minmod_dx(gp1, g, gm1);
        WmhR = g - 0.5*minmod->minmod_dx(gp1, g, gm1);
        WmhL = gm1 + 0.5*minmod->minmod_dx(g, gm1, gm2);

        //a = fabs(grid_pt->u[rk_flag][direc]);
        a = fabs(vis_array[idx][16])/vis_array[idx][15];
        am1 = (fabs(vis_nbr_x[1][16])/vis_nbr_x[1][15]);
        ap1 = (fabs(vis_nbr_x[2][16])/vis_nbr_x[2][15]);

        ax = maxi(a, ap1);
        HWph = ((uWphR + uWphL) - ax*(WphR - WphL))*0.5;

        ax = maxi(a, am1);
        HWmh = ((uWmhR + uWmhL) - ax*(WmhR - WmhL))*0.5;
        
        HW = (HWph - HWmh)/DATA_ptr->delta_x/taufactor;
            
        // make partial_i (u^i Wmn)
        sum += -HW;
        
        // y-direction
        taufactor = 1.0;
        /* Get_uWmns */
        g = vis_array[idx][idx_1d];
        f = g*vis_array[idx][17];
        g *= vis_array[idx][15];
           
        gp2 = vis_nbr_y[3][idx_1d];
        fp2 = gp2*vis_nbr_y[3][17];
        gp2 *= vis_nbr_y[3][15];
        
        gp1 = vis_nbr_y[2][idx_1d];
        fp1 = gp1*vis_nbr_y[2][17];
        gp1 *= vis_nbr_y[2][15];
        
        gm1 = vis_nbr_y[1][idx_1d];
        fm1 = gm1*vis_nbr_y[1][17];
        gm1 *= vis_nbr_y[1][15];
        
        gm2 = vis_nbr_y[0][idx_1d];
        fm2 = gm2*vis_nbr_y[0][17];
        gm2 *= vis_nbr_y[0][15];

        // MakeuWmnHalfs uWmn
        uWphR = fp1 - 0.5*minmod->minmod_dx(fp2, fp1, f);
        uWphL = f + 0.5*minmod->minmod_dx(fp1, f, fm1);
        uWmhR = f - 0.5*minmod->minmod_dx(fp1, f, fm1);
        uWmhL = fm1 + 0.5*minmod->minmod_dx(f, fm1, fm2);

        // just Wmn
        WphR = gp1 - 0.5*minmod->minmod_dx(gp2, gp1, g);
        WphL = g + 0.5*minmod->minmod_dx(gp1, g, gm1);
        WmhR = g - 0.5*minmod->minmod_dx(gp1, g, gm1);
        WmhL = gm1 + 0.5*minmod->minmod_dx(g, gm1, gm2);

        //a = fabs(grid_pt->u[rk_flag][direc]);
        a = fabs(vis_array[idx][17])/vis_array[idx][15];
        am1 = (fabs(vis_nbr_y[1][17])/vis_nbr_y[1][15]);
        ap1 = (fabs(vis_nbr_y[2][17])/vis_nbr_y[2][15]);
        ax = maxi(a, ap1);
        HWph = ((uWphR + uWphL) - ax*(WphR - WphL))*0.5;
        ax = maxi(a, am1);
        HWmh = ((uWmhR + uWmhL) - ax*(WmhR - WmhL))*0.5;
        HW = (HWph - HWmh)/DATA_ptr->delta_x/taufactor;
            
        // make partial_i (u^i Wmn)
        sum += -HW;
        
        // eta-direction
        taufactor = tau;
        /* Get_uWmns */
        g = vis_array[idx][idx_1d];
        f = g*vis_array[idx][18];
        g *= vis_array[idx][15];
           
        gp2 = vis_nbr_eta[3][idx_1d];
        fp2 = gp2*vis_nbr_eta[3][18];
        gp2 *= vis_nbr_eta[3][15];
        
        gp1 = vis_nbr_eta[2][idx_1d];
        fp1 = gp1*vis_nbr_eta[2][18];
        gp1 *= vis_nbr_eta[2][15];
        
        gm1 = vis_nbr_eta[1][idx_1d];
        fm1 = gm1*vis_nbr_eta[1][18];
        gm1 *= vis_nbr_eta[1][15];
        
        gm2 = vis_nbr_eta[0][idx_1d];
        fm2 = gm2*vis_nbr_eta[0][18];
        gm2 *= vis_nbr_eta[0][15];

        // MakeuWmnHalfs uWmn
        uWphR = fp1 - 0.5*minmod->minmod_dx(fp2, fp1, f);
        uWphL = f + 0.5*minmod->minmod_dx(fp1, f, fm1);
        uWmhR = f - 0.5*minmod->minmod_dx(fp1, f, fm1);
        uWmhL = fm1 + 0.5*minmod->minmod_dx(f, fm1, fm2);

        // just Wmn
        WphR = gp1 - 0.5*minmod->minmod_dx(gp2, gp1, g);
        WphL = g + 0.5*minmod->minmod_dx(gp1, g, gm1);
        WmhR = g - 0.5*minmod->minmod_dx(gp1, g, gm1);
        WmhL = gm1 + 0.5*minmod->minmod_dx(g, gm1, gm2);

        //a = fabs(grid_pt->u[rk_flag][direc]);
        a = fabs(vis_array[idx][18])/vis_array[idx][15];
        am1 = (fabs(vis_nbr_eta[1][17])/vis_nbr_eta[1][18]);
        ap1 = (fabs(vis_nbr_eta[2][18])/vis_nbr_eta[2][18]);
        ax = maxi(a, ap1);
        HWph = ((uWphR + uWphL) - ax*(WphR - WphL))*0.5;
        ax = maxi(a, am1);
        HWmh = ((uWmhR + uWmhL) - ax*(WmhR - WmhL))*0.5;
        HW = (HWph - HWmh)/DATA_ptr->delta_x/taufactor;
            
        // make partial_i (u^i Wmn)
        sum += -HW;
    
        w_rhs[mu][nu] = sum*(DATA_ptr->delta_tau);
    }
    return(1);
}

double Diss::get_temperature_dependent_eta_s(double T) {
    double Ttr = 0.18/hbarc;  // phase transition temperature
    double Tfrac = T/Ttr;
    double shear_to_s;
    if (T < Ttr) {
        shear_to_s = (DATA_ptr->shear_to_s + 0.0594*(1. - Tfrac)
                      + 0.544*(1. - Tfrac*Tfrac));
    } else {
        shear_to_s = (DATA_ptr->shear_to_s + 0.288*(Tfrac - 1.) 
                      + 0.0818*(Tfrac*Tfrac - 1.));
    }
    return(shear_to_s);
}

double Diss::get_temperature_dependent_zeta_s(double temperature) {
    // T dependent bulk viscosity from Gabriel
    /////////////////////////////////////////////
    //           Parametrization 1             //
    /////////////////////////////////////////////
    double Ttr=0.18/0.1973;
    double dummy=temperature/Ttr;
    double A1=-13.77, A2=27.55, A3=13.45;
    double lambda1=0.9, lambda2=0.25, lambda3=0.9, lambda4=0.22;
    double sigma1=0.025, sigma2=0.13, sigma3=0.0025, sigma4=0.022;
 
    double bulk = A1*dummy*dummy + A2*dummy - A3;
    if (temperature < 0.995*Ttr) {
        bulk = (lambda3*exp((dummy-1)/sigma3)
                + lambda4*exp((dummy-1)/sigma4) + 0.03);
    }
    if (temperature > 1.05*Ttr) {
        bulk = (lambda1*exp(-(dummy-1)/sigma1)
                + lambda2*exp(-(dummy-1)/sigma2) + 0.001);
    }
    
    /////////////////////////////////////////////
    //           Parametrization 2             //
    /////////////////////////////////////////////
    //double Ttr=0.18/0.1973;
    //double dummy=temperature/Ttr;
    //double A1=-79.53, A2=159.067, A3=79.04;
    //double lambda1=0.9, lambda2=0.25, lambda3=0.9, lambda4=0.22;
    //double sigma1=0.025, sigma2=0.13, sigma3=0.0025, sigma4=0.022;

    //bulk = A1*dummy*dummy + A2*dummy - A3;

    //if (temperature < 0.997*Ttr) {
    //    bulk = (lambda3*exp((dummy-1)/sigma3)
    //            + lambda4*exp((dummy-1)/sigma4) + 0.03);
    //}
    //if (temperature > 1.04*Ttr) {
    //    bulk = (lambda1*exp(-(dummy-1)/sigma1)
    //            + lambda2*exp(-(dummy-1)/sigma2) + 0.001);
    //}

    ////////////////////////////////////////////
    //           Parametrization 3            //
    ////////////////////////////////////////////
    //double Ttr=0.18/0.1973;
    //double dummy=temperature/Ttr;
    //double lambda1=0.9, lambda2=0.25, lambda3=0.9, lambda4=0.22;
    //double sigma1=0.025, sigma2=0.13, sigma3=0.0025, sigma4=0.022;
    
    //if (temperature<0.99945*Ttr) {
    //    bulk = (lambda3*exp((dummy-1)/sigma3)
    //            + lambda4*exp((dummy-1)/sigma4) + 0.03);
    //}
    //if (temperature>0.99945*Ttr) {
    //    bulk = 0.901*exp(14.5*(1.0-dummy)) + 0.061/dummy/dummy;
    //}

    return(bulk);
}

