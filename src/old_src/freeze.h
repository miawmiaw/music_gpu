// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_FREEZE_H_
#define SRC_FREEZE_H_

#define NUMDECAY 2000
#define MAXINTV 2000000     // size of arry for Montecarlo numbers
                            // note that I modified pdg05.dat
                            // by replacing 90...... 9......*/
#define MHALF (MAXINTV/2)
#define NY 200              /* size of arry for storage of the y-spectrum */
#define NPT 100             /* size of arry for storage of the pt-spectrum */
#define NPHI 100            /* size of arry for storage of the phi-spectrum */
#define NPHI1 NPHI + 1
#define PTS3 12             /* normalization of 3-body decay */
#define PTS4 12             /* inv. mass integral 3-body    */

#define PTN1 8
#define PTN2 8              /* 2-body between the poles */

#define PTCHANGE    1.

#include <unistd.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <string>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "./data.h"
#include "./util.h"
#include "./eos.h"

const int nharmonics = 8;   // calculate up to maximum harmonic (n-1)
                            // -- for nharmonics = 8, calculate from v_0 o v_7
const int etasize = 100;    // max number of points in eta for array
const int ptsize = 100;     // max number of points in pt for array

class Freeze {
 private:
    typedef struct particle {
        int number;
        char *name;
        double mass;
        double width;
        int degeneracy;
        int baryon;
        int strange;
        int charm;
        int bottom;
        int isospin;
        double charge;
        int decays;
        int stable;
        int ny;
        int npt;
        int nphi;
        double phimin;
        double phimax;
        double ymax;
        double deltaY;
        double  resCont[NY][NPT][NPHI];
        double  dNdydptdphi[NY][NPT][NPHI+1];
        double  dNdydptdphi2[NY][NPT][NPHI+1];
        double  mt[NPT];    /* mt values for spectrum */
        double  pt[NPT];    /* pt values for spectrum */
        double  y[NY];      /* y values for spectrum */
        double  slope;      /* assymtotic slope of pt-spectrum */
        double muAtFreezeOut;   // the chemical potential at freeze-out
                                // for the partial chemical equilibrium
    } Particle;

    struct de {
        int  reso;       /* Montecarlo number of decaying resonance */
        int  numpart;    /* number of daughter particles after decay */
        double branch;     /* branching ratio */
        int    part[5];    /* array of daughter particles Montecarlo number */
    } decay[NUMDECAY];

    typedef struct pblockN {
        double pt, mt, y, e, pl;    /* pt, mt, y of decay product 1 */
        double phi;
        double m1, m2, m3;          /* masses of decay products     */
        double mr;                  /* mass of resonance            */
        double costh, sinth;
        double e0, p0;
        int res_num;                /* Montecarlo number of the Res. */
    } pblockN;

    typedef struct nblock {
        double a, b, c, d;
    } nblock;         /* for normalisation integral of 3-body decays */

    typedef struct surfaceElement {
        double x[4];        // position in (tau, x, y, eta)
        double cosh_eta_s, sinh_eta_s;  // caching the sinh and cosh of eta_s
        double s[4];        // hypersurface vector in (tau, x, y, eta)
        double u[4];        // flow velocity in (tau, x, y, eta)
        double W[4][4];     // W^{\mu\nu}
        double q[4];        // baryon diffusion current
        double pi_b;        // bulk pressure
        double rho_B;       // net baryon density
        double epsilon_f;
        double T_f;
        double mu_B;
        double eps_plus_p_over_T_FO;  // (energy_density+pressure)/temperature
    } SurfaceElement;
    SurfaceElement *surface;
    Particle *particleList;
    int NCells;
    int decayMax, particleMax;

    double ***sumYPtPhi;
    int *partid;
    // array for converting Montecarlo numbers in internal numbering of
    // the resonances
    double *phiArray;
    // Int *integral;
    Util *util;
    int pseudofreeze;

    InitData* DATA_ptr;
    int *charged_hadron_list;
    int charged_hadron_list_length;

    int deltaf_qmu_coeff_table_length_T;
    int deltaf_qmu_coeff_table_length_mu;
    double delta_qmu_coeff_table_T0, delta_qmu_coeff_table_mu0;
    double delta_qmu_coeff_table_dT, delta_qmu_coeff_table_dmu;
    double **deltaf_qmu_coeff_tb;

    int deltaf_coeff_table_14mom_length_T;
    int deltaf_coeff_table_14mom_length_mu;
    double delta_coeff_table_14mom_T0, delta_coeff_table_14mom_mu0;
    double delta_coeff_table_14mom_dT, delta_coeff_table_14mom_dmu;
    double **deltaf_coeff_tb_14mom_DPi, **deltaf_coeff_tb_14mom_BPi;
    double **deltaf_coeff_tb_14mom_BPitilde;
    double **deltaf_coeff_tb_14mom_BV, **deltaf_coeff_tb_14mom_DV;
    double **deltaf_coeff_tb_14mom_Bpi_shear;

    int bulk_deltaf_kind;

 public:
    Freeze(InitData* DATA_in);  // constructor
    ~Freeze();  // destructor

    double gauss(int n, double (Freeze::*f)(double, void *),
                 double xlo, double xhi, void *optvec);
    double riemannsum(int n, double (Freeze::*f)(double, void *),
                      double xlo, double xhi, void *optvec);
    void ReadParticleData(InitData *DATA, EOS *eos);
    void ReadFreezeOutSurface(InitData *DATA);
    void ReadSpectra(InitData *DATA);
    void Read3Spectra(InitData *DATA);
    void ReadSingleSpectrum(InitData* DATA);
    void ReadFullSpectra(InitData *DATA);
    void ReadFullSpectra2(InitData *DATA);
    void ComputeAveragePT(int number, double ptmax);
    void ComputeChargedHadrons(InitData* DATA, double ptmax);
    void Compute3ChargedHadrons(InitData *DATA, double ptmax);
    void ComputeCorrelations(InitData* DATA, double ptmax);
    double summation(double px, double py, double y, double m, int deg,
                     int baryon, double mu, InitData *DATA);
    double summation3(double px, double py, double y, double m, int deg,
                      int baryon, double mu, InitData *DATA);
    void ComputeParticleSpectrum(InitData *DATA, int number, double ptmax,
                                 int anti, int iptmax, int iphimax, int size,
                                 int rank);
    void OutputFullParticleSpectrum(InitData *DATA, int number, double ptmax,
                                    int anti, int full);
    void load_deltaf_qmu_coeff_table(string filename);
    void load_deltaf_qmu_coeff_table_14mom(string filename);
    void getbulkvisCoefficients(double Tdec, double* bulkvisCoefficients);

    // -------------------------------------------------------------------------
    // the following routines are adapted from the public version of
    /* ... the resonance decay calculation using the output  */
    /* generated by the hydrodynamical code azhydro0p2. 
    /* The majority of the code was  */
    /* developed by Josef Sollfrank and Peter Kolb.
    /* Additional work and comments  */
    /* were added by Evan Frodermann, September 2005. */
    /* Please refer to the papers  */
    /* J. Sollfrank, P. Koch, and U. Heinz, Phys. Lett B 252 (1990) and  */
    /* J. Sollfrank, P. Koch, and U. Heinz, Z. Phys. C 52 (1991)  */
    /* for a description of the formalism utilized in this program. */

    // this computes "Q(m_R,m_1,m_2,m_3)"
    double norm3int(double x, void *paranorm);

    double Edndp3(double yr, double ptr, double phirin, int res_num);
    double dnpir2N(double phi, void *para1);
    double dnpir1N(double costh, void* para1);
    double dn2ptN(double w2, void* para1);
    double dn3ptN(double x, void* para1);
    double Edndp3_2bodyN(double y, double pt, double phi, double m1, double m2,
                         double mr, int res_num);
    double Edndp3_3bodyN(double y, double pt, double phi, double m1, double m2,
                         double m3, double mr, double norm3, int res_num);
    void add_reso(int pn, int pnR, int k, int j);
    void cal_reso_decays(int maxpart, int maxdecay, int bound, int mode);
    // -------------------------------------------------------------------------
    int countLines(std::istream& in);
    void checkForReadError(FILE *file, const char* name);
    void CooperFrye(int particleSpectrumNumber, int mode, InitData *DATA,
                    EOS *eos, int size, int rank);

    // When the pseudorapidity mode is used
    void ReadFSpectra_pseudo(InitData *DATA);
    void ReadSpectra_pseudo(InitData* DATA, int full, int verbose);
    void ComputeParticleSpectrum_pseudo(InitData *DATA, int number, int anti, int size, int rank);
  void ComputeParticleSpectrum_pseudo_improved(InitData *DATA, int number, int size, int rank);
  void OutputFullParticleSpectrum_pseudo(InitData *DATA, int number, int anti, int full);
  void CooperFrye_pseudo(int particleSpectrumNumber, int mode, InitData *DATA, EOS *eos, int size, int rank);
  double    Edndp3_pseudo(double yr, double ptr, double phirin, int res_num);
  double dnpir2N_pseudo (double phi, void *para1);    
  double dnpir1N_pseudo (double costh, void* para1);
  double dn2ptN_pseudo (double w2, void* para1);
  double dn3ptN_pseudo (double x, void* para1);
  double Edndp3_2bodyN_pseudo (double y, double pt, double phi, double m1, double m2, double mr, int res_num);
  double Edndp3_3bodyN_pseudo (double y, double pt, double phi, double m1, double m2,
            double m3, double mr, double norm3, int res_num);
  void add_reso_pseudo (int pn, int pnR, int k, int j, int pseudofreeze);
  void cal_reso_decays_pseudo (int maxpart, int maxdecay, int bound, int mode, int pseudofreeze);
  double Rap(double eta, double pt, double m);
  double PseudoRap(double y, double pt, double m);
  double dydeta(double eta, double pt, double m);
//   void pt_integrated_flow(InitData *DATA, int number, double minpt, double maxpt, double ****vn);
  void pt_integrated_flow(InitData *DATA, int number, double minpt, double maxpt, double vn[nharmonics][2][etasize]);
  void eta_integrated_flow(InitData *DATA, int number, double mineta, double maxeta, double vn[nharmonics][2][ptsize]);
  void y_integrated_flow(InitData *DATA, int number, double mineta, double maxeta, double vn[nharmonics][2][ptsize]);
  void rapidity_integrated_flow(InitData *DATA, int number, int yflag, double minrap, double maxrap, double vn[nharmonics][2][ptsize]);
  void pt_and_eta_integrated_flow(InitData *DATA, int number, double minpt, double maxpt, double mineta, double maxeta, double vn[8][2]);
  void pt_and_eta_integrated_flow2(InitData *DATA, int number, double minpt, double maxpt, double mineta, double maxeta, double vn[8][2]);
  void pt_and_y_integrated_flow(InitData *DATA, int number, double minpt, double maxpt, double miny, double maxy, double vn[8][2]);
  void pt_and_rapidity_integrated_flow(InitData *DATA, int number, double minpt, double maxpt, int yflag, double minrap, double maxrap, double vn[8][2]);
  void weighted_v1(InitData *DATA, int number, double minpt, double maxpt, int yflag, double minrap, double maxrap, double vn[2], int ch);
  double get_yield(InitData *DATA, int number, double minpt, double maxpt, double miny, double maxy);
  double get_N(InitData *DATA, int number, double minpt, double maxpt, int yflag, double minrap, double maxrap);
  double get_meanpt(InitData *DATA, int number, double minpt, double maxpt, int yflag, double minrap, double maxrap);
  double get_vn(InitData *DATA, int number, double minpt, double maxpt, int yflag, double minrap, double maxrap, int n);
  double get_psi_n(InitData *DATA, int number, double minpt, double maxpt, int yflag, double minrap, double maxrap, int n);
  double get_Nch(InitData *DATA, double minpt, double maxpt, int yflag, double minrap, double maxrap);
  double get_vn_ch(InitData *DATA, double minpt, double maxpt, int yflag, double minrap, double maxrap, int n, double* vn_results);
  double get_psi_n_ch(InitData *DATA, double minpt, double maxpt, int yflag, double minrap, double maxrap, int n);
  double get_weighted_v1(InitData *DATA, int number, double minpt, double maxpt, int yflag, double minrap, double maxrap, int ch);
  double get_weighted_psi1(InitData *DATA, int number, double minpt, double maxpt, int yflag, double minrap, double maxrap, int ch);

  void OutputDifferentialFlowAtMidrapidity(InitData *DATA, int number, int full);
  void OutputDifferentialFlowNearMidrapidity(InitData *DATA, int number, int full);
  void OutputIntegratedFlow_vs_y(InitData *DATA, int number, int full, double pT_min, double pT_max);
  void Output_charged_IntegratedFlow(InitData *DATA, double pT_min, double pT_max, double eta_min, double eta_max);
  void Output_charged_hadrons_eta_differential_spectra(InitData *DATA, int full, double pT_min ,double pT_max);
  void Output_charged_hadrons_pT_differential_spectra(InitData *DATA, int full, double eta_min, double eta_max);
  void Output_midrapidity_hadrons_spectra(InitData *DATA, int full, const int * hadron_list, int nb_hadrons);

   double get_deltaf_qmu_coeff(double T, double muB);
   double get_deltaf_coeff_14moments(double T, double muB, double type);
};
#endif
  
