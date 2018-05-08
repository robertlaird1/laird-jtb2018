// sequential_game
// 
// Runs the sequential, one-shot Prisoner's Dilemma and Snowdrift games (spatial or aspatial)
// Matlab call: "[tsabund, tsoutcome, snapshot] = sequential_game(R, S, T, P, L, startcond, beta, mu, gmax, seq, sp, snapon, snapint)" 
// Robert Laird, University of Lethbridge, Sept. 2, 2015 (robert.laird@uleth.ca)
// Includes original 'Mersenne Twister' code by Takuji Nishimura and Makoto Matsumoto (2002); included and modified as allowed in their code (see Appendix 5, below)

#include <mex.h>
#include <math.h>

using namespace std;

// 1. SUB-ROUTINE AND FUNCTION PROTOTYPES

double twist(void); // function: returns a random double on [0, 1) using Mersenne Twister

// 2. MAIN FUNCTION

void sim(
    // output (lhs) variables
    double tsabund_lhs[],       // time series: abundance values
    double tsoutcome_lhs[],     // time series: outcomes of interactions
    double snapshot_lhs[],      // snapshot of the population

    // input (rhs) variables
    double R_rhs[],             // game characteristics: reward payoff
    double S_rhs[],             // game characteristics: sucker's payoff
    double T_rhs[],             // game characteristics: temptation payoff
    double P_rhs[],             // game characteristics: punishment payoff 
    double L_rhs[],             // population structure: square-root of population size
    double startcond_rhs[],     // simulation details: starting conditions    
    double beta_rhs[],          // simulation details: beta parameter (shape parameter; applies to Fermi rule only; set beta = -1 to invoke replicator rule)
    double mu_rhs[],            // simulation details: mutation rate
    double gmax_rhs[],          // simulation details: number of model generations
    double seq_rhs[],           // simulation details: sequential (1) or simultaneous (0)
    double sp_rhs[],            // simulation details: spatial (1) or well-mixed (0)    
    double snapon_rhs[],        // output options: record snapshots (on/off switch)
    double snapint_rhs[]        // output options: snapshot interval (in generations)
    )        
{

    // 2.1. DECLARE / INITIALIZE VARIABLES DIRECTLY RELATED TO RHS ARGUMENTS (PLUS THE ALL-PURPOSE COUNTER, ii)
    
    const double R = R_rhs[0]; // reward payoff for mutual cooperation
    const double S = S_rhs[0]; // sucker's payoff for unilateral cooperation with a defector
    const double T = T_rhs[0]; // temptation payoff for unilateral defection against a cooperator
    const double P = P_rhs[0]; // punishment payoff for mutual defection
    const int L = static_cast<int>(L_rhs[0]); // square-root of the population size
    double startcond[8]; // starting conditions: relative requence of the strategies
    for (int j = 0; j < 8; j++) {startcond[j] = startcond_rhs[j];}
    const double beta = beta_rhs[0]; // beta parameter (shape parameter; applies to Fermi rule only; set to -1 to invoke replicator rule)
    const double mu = mu_rhs[0]; // mutation rate
    const int gmax = static_cast<int>(gmax_rhs[0]); // maximum number of generations
    const bool seq = (seq_rhs[0] != 0); // sequential game, off (0) or on (1)
    const bool sp = (sp_rhs[0] != 0); // spatial arena, off (0) or on (1)
    const bool snapon = (snapon_rhs[0] != 0); // movie frame data, off (0) or on (1)
    const int snapint = static_cast<int>(snapint_rhs[0]); // move frame interval, in generations
    int ii = 0; // an all-purpose counter; used multiple times in this code
        
    // 2.2. DECLARE / INITIALIZE VARIABLES INDIRECTLY RELATED TO RHS ARGUMENTS (I.E., SET UP INITIAL CONDITIONS AND DATA RECORDERS)
    
    // 2.2.1. startcondPMF, starting conditions scaled so that the elements sum to 1 (PMF = probability mass function)
    // 2.2.2. startcondCDF, sum of all the elements in startcondPMF up to and including the equivalent value of the index (CDF = cumulative density function)
    if (!seq) { // in the non-sequential version, calls to include strategies 1-6 are ignored; only 0 and 7 are possible
        for (int j = 1; j < 7; j++) {startcond[j] = 0;} // set starting abund off invalid strategies to 0
    }
    double startcondPMF[8];
    double startcondCDF[8];
    double startcondsum = 0.0; // sum of elements in startcond
    for (int j = 0; j < 8; j++) {startcondsum += startcond[j];}
    double startcondPMFsum = 0.0; // sum of elements in startcondPMF
    for (int j = 0; j < 8; j++) {
        if (startcondsum == 0) { // i.e., if the user put 0s for every strategy in startcond
            if (seq) { // sequential version of the model
                startcondPMF[j] = 0.125;
            }
            else { // non-sequential version of the model
                if ((j == 0) || (j == 7)) { // only two valid strategies in this version
                    startcondPMF[j] = 0.5; 
                }
                else {
                    startcondPMF[j] = 0;
                }
            }
        }
        else {
            startcondPMF[j] = startcond[j]/startcondsum; // proportions for every strategy at the start
        }
        startcondPMFsum += startcondPMF[j];
        startcondCDF[j] = startcondPMFsum;
    }
    
    // 2.2.3. world, the competive arena (square L x L lattice, with toroidal boundaries, of degree 4 (von Neumann neighbourhood))
    // 2.2.4. abund, the abundance of the eight strategies
    int* world = new int[L*L]; // L rows by L columns
    int abund[8]; // running abundance of the eight strategies (0 to 7)
    for (int j = 0; j < 8; j++) {abund[j] = 0;} // initialize abundance
    double rd; // random double on the interval [0, 1)
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            
            // 2.2.4.1. randomly choose a strategy, scaled by the starting conditions
            rd = twist(); // random double on the interval [0, 1)
            for (int k = 0; k < 8; k++) {
                if (k == 0) {
                    if (rd < startcondCDF[k]) {
                        world[i*L + j] = k; // put an individual with strategy k (k = 0 in this case) in the world
                        abund[k]++; // increment the abundance of strategy k (k = 0 in this case)
                    }
                }
                else {
                    if ((rd >= startcondCDF[k - 1]) && (rd < startcondCDF[k])) {
                        world[i*L + j] = k; // put an individual with strategy k in the world
                        abund[k]++; // increment the abundance of strategy k
                    }                
                }            
            }
            
        }
    }
    
    // 2.2.5. payoff1st, the payoff matrix for an individual who plays first
    const double payoff1st[64] = {R, R, S, S, R, R, S, S,
                                  R, R, S, S, R, R, S, S,
                                  R, R, S, S, R, R, S, S,
                                  R, R, S, S, R, R, S, S,
                                  T, P, T, P, T, P, T, P,
                                  T, P, T, P, T, P, T, P,
                                  T, P, T, P, T, P, T, P,
                                  T, P, T, P, T, P, T, P};
                               
    // 2.2.6. payoff2nd, the payoff matrix for an individual who plays second
    const double payoff2nd[64] = {R, R, R, R, S, S, S, S,
                                  R, R, R, R, P, P, P, P,
                                  T, T, T, T, S, S, S, S,
                                  T, T, T, T, P, P, P, P,
                                  R, R, R, R, S, S, S, S,
                                  R, R, R, R, P, P, P, P,
                                  T, T, T, T, S, S, S, S,
                                  T, T, T, T, P, P, P, P};
                                  
    // 2.2.7. alpha, the difference between the highest and lowest payoff                              
    double alpha; // difference between highest and lowest payoff
    double maxpay = R; // maximum payoff (initialized as R for comparison purposes)
    double minpay = R; // minimum payoff (initialized as R for comparison purposes)
    if (S > maxpay) {maxpay = S;}
    if (T > maxpay) {maxpay = T;} // i.e., as in the Prisoner's Dilemma and Snowdrift games
    if (P > maxpay) {maxpay = P;}
    if (S < minpay) {minpay = S;} // i.e., as in the Prisoner's Dilemma
    if (T < minpay) {minpay = T;}
    if (P < minpay) {minpay = P;} // i.e., as in the Snowdrift game
    alpha = maxpay - minpay;    
                                  
    // 2.2.8. tsabund, a data recorder for the time series of abundance
    int* tsabund = new int[(gmax + 1)*8]; // (gmax + 1) rows by 8 columns
    for (int j = 0; j < 8; j++) {
        tsabund[0*8 + j] = abund[j]; // record initial abundances in row 0
    }
    
    // 2.2.9. tsoutcome, a data recorder for the time series of outcomes of interactions
    int* tsoutcome = new int[gmax*4]; // gmax rows by 4 columns
    
    // 2.2.10. define and initialize snapshot
    int nsnap = 1; // number of snapshots to record, if recorded (otherwise, dummy entry of 1)
    int snapsize = 1; // size of the entire set of snapshots, if recorded (otherwise, dummy entry of 1)
    if (snapon) { // record snapshots
        nsnap = static_cast<int>(gmax/snapint + 1); // number of snapshots to record (will truncate/round down to int)
        snapsize = nsnap*L*L; // size of the entire set of snapshots
    }
    int* snapshot = new int[snapsize]; // data for creating snapshots for (tag)stratagies (nsnap rows by L*L columns)
    
    // 2.2.11. record first snapshot
    int snapID = 0; // gives the snapID number for the set of snapshots
    if (snapon) { // i.e., if recording a snapshot  
        for (int j = 0; j < L*L; j++) {
            snapshot[snapID*L*L + j] = world[j]; // add row 0 (snapshot of initial conditions)
        }
    }
        
    // 2.3 SIMULATION
    int outcome[4]; // records outcome of every interaction (position 0 records CC interactions, 1 records CD, 2 records DC, 3 records DD)
    int fr; // focal row (integer between 0 and L - 1)
    int fc; // focal col (integer between 0 and L - 1)
    int cd; // competitor direction (random integer between 0 and 3)
    int cr; // comp row (integer between 0 and L - 1)
    int cc; // comp col (integer between 0 and L - 1)
    double fpay; // focal's average payoff
    double cpay; // comps's average payoff
    int fnr; // focal neighbour row (integer between 0 and L - 1)
    int fnc; // focal neighbour col (integer between 0 and L - 1)
    int cnr; // comp neighbour row (integer between 0 and L - 1)
    int cnc; // comp neighbour col (integer between 0 and L - 1)
    int mutant; // strategy ID of mutant
    for (int g = 0; g < gmax; g++) { // generation loop
        
        for (int j = 0; j < 4; j++) {outcome[j] = 0;} // reset the outcome counter
        
        for (int t = 0; t < L*L; t++) { // time loop
        
            // 2.3.1. determine fr and fc, the focal individual's row and column
            fr = static_cast<int>(twist()*L); // focal row (random integer between 0 and L - 1)
            fc = static_cast<int>(twist()*L); // focal col (random integer between 0 and L - 1)
            
            // 2.3.2. determine cd, the competitor's direction, relative to focal
            if (sp) {cd = static_cast<int>(twist()*4);} // competitor direction for spatial version (random integer between 0 and 3)
            
            // 2.3.3. determine cr and cc, the competitor individual's row and column (accounting for wrapping boundaries)
            if (sp) { // spatial version of the model
                if (cd == 0) { // comp is North of focal
                    if (fr == 0) {cr = L - 1;} else {cr = fr - 1;} // comp row (integer between 0 and L - 1)
                    cc = fc; // comp col (integer between 0 and L - 1)
                }
                else if (cd == 1) { // comp is West of focal
                    cr = fr; // comp row (integer between 0 and L - 1)
                    if (fc == 0) {cc = L - 1;} else {cc = fc - 1;} // comp col (integer between 0 and L - 1)
                }
                else if (cd == 2) { // comp is East of focal
                    cr = fr; // comp row (integer between 0 and L - 1)
                    if (fc == L - 1) {cc = 0;} else {cc = fc + 1;} // comp col (integer between 0 and L - 1)
                }
                else { // comp is South of focal
                    if (fr == L - 1) {cr = 0;} else {cr = fr + 1;} // comp row (integer between 0 and L - 1)
                    cc = fc; // comp col (integer between 0 and L - 1)
                }
            }
            else { // aspatial version of the model
                do {
                    cr = static_cast<int>(twist()*L); // comp row (random integer between 0 and L - 1)
                    cc = static_cast<int>(twist()*L); // comp col (random integer between 0 and L - 1) 
                } while ((fr == cr) && (fc == cc)); // to make sure comp is different than focal
            }
            
            // 2.3.4. determine payoffs when focal and comp play their four nearest neighbours (spatial) or four random individuals (aspatial); determine outcome of interactions
            fpay = 0.0; // initialize focal's payoff
            cpay = 0.0; // initialize comp's payoff
            for (int nd = 0; nd < 4; nd++) { // directions of four nearest neighbours (or four random individuals, as the case may be)
                if (sp) { // spatial version of the model
                    if (nd == 0) { // North neighbours
                        if (fr == 0) {fnr = L - 1;} else {fnr = fr - 1;} // focal neighbour row (integer between 0 and 1)
                        fnc = fc; // focal neighbour col (integer between 0 and 1)
                        if (cr == 0) {cnr = L - 1;} else {cnr = cr - 1;} // comp neighbour row (integer between 0 and 1)
                        cnc = cc; // comp neighbour col (integer between 0 and 1)
                    }
                    else if (nd == 1) { // West neighbours
                        fnr = fr; // focal neighbour row (integer between 0 and 1)
                        if (fc == 0) {fnc = L - 1;} else {fnc = fc - 1;} // focal neighbour col (integer between 0 and 1)
                        cnr = cr; // comp neighbour row (integer between 0 and 1)
                        if (cc == 0) {cnc = L - 1;} else {cnc = cc - 1;} // comp neighbour col (integer between 0 and 1)
                    }
                    else if (nd == 2) { // East neighbours
                        fnr = fr; // focal neighbour row (integer between 0 and 1)
                        if (fc == L - 1) {fnc = 0;} else {fnc = fc + 1;}  // focal neighbour col (integer between 0 and 1)
                        cnr = cr; // comp neighbour row (integer between 0 and 1)
                        if (cc == L - 1) {cnc = 0;} else {cnc = cc + 1;} // comp neighbour col (integer between 0 and 1)
                    }
                    else { // South neighbours
                        if (fr == L - 1) {fnr = 0;} else {fnr = fr + 1;} // focal neighbour row (integer between 0 and 1)
                        fnc = fc; // focal neighbour col (integer between 0 and 1)
                        if (cr == L - 1) {cnr = 0;} else {cnr = cr + 1;} // comp neighbour row (integer between 0 and 1)
                        cnc = cc; // comp neighbour col (integer between 0 and 1)
                    }
                }
                else { // aspatial version of the model
                    do {
                        fnr = static_cast<int>(twist()*L); // focal neighbour row (random integer between 0 and L - 1)
                        fnc = static_cast<int>(twist()*L); // focal neighbour col (random integer between 0 and L - 1)
                    } while ((fr == fnr) && (fc == fnc)); // to make sure focal neighbour is different than focal
                    do {
                        cnr = static_cast<int>(twist()*L); // comp neighbour row (random integer between 0 and L - 1)
                        cnc = static_cast<int>(twist()*L); // comp neighbour col (random integer between 0 and L - 1)
                    } while ((cr == cnr) && (cc == cnc)); // to make sure comp neighbour is different than comp
                }
                // note that order doesn't matter in the non-seq version when only strats 0 and 7 are possible
                if (twist() < 0.5) { // focal goes first against its neighbour (order determined randomly)
                    fpay += payoff1st[world[fr*L + fc]*8 + world[fnr*L + fnc]]; // add to focal's payoff
                    if (world[fr*L + fc] < 4) { // focal cooperates first 
                        if ((world[fnr*L + fnc] == 0) || (world[fnr*L + fnc] == 1) || (world[fnr*L + fnc] == 4) || (world[fnr*L + fnc] == 5)) { // focal neighbour responds with cooperation
                            outcome[0]++; // increment outcome 0 (CC)
                        }
                        else { // focal neighbour responds with defection 
                            outcome[1]++; // increment outcome 1 (CD)
                        }
                    }
                    else { // focal defects first
                        if ((world[fnr*L + fnc] == 0) || (world[fnr*L + fnc] == 2) || (world[fnr*L + fnc] == 4) || (world[fnr*L + fnc] == 6)) { // focal neighbour responds with cooperation
                            outcome[2]++; // increment outcome 2 (DC) 
                        }
                        else { // focal neighbour responds with defection
                            outcome[3]++; // increment outcome 3 (DD)
                        }
                    }
                }
                else { // focal goes second against its neighbour (order determined randomly)
                    fpay += payoff2nd[world[fr*L + fc]*8 + world[fnr*L + fnc]]; // add to focal's payoff
                    if (world[fnr*L + fnc] < 4) { // focal neighbour cooperates first
                        if ((world[fr*L + fc] == 0) || (world[fr*L + fc] == 1) || (world[fr*L + fc] == 4) || (world[fr*L + fc] == 5)) { // focal responds with cooperation
                            outcome[0]++; // increment outcome 0 (CC)
                        }
                        else { // focal responds with defection
                            outcome[1]++; // increment outcome 1 (CD)
                        }
                    }
                    else { // focal neighbour defects first
                        if ((world[fr*L + fc] == 0) || (world[fr*L + fc] == 2) || (world[fr*L + fc] == 4) || (world[fr*L + fc] == 6)) { // focal responds with cooperation
                            outcome[2]++; // increment outcome 2 (DC)
                        }
                        else { // focal responds with defection
                            outcome[3]++; // increment outcome 3 (DD)
                        }
                    }  
                }
                if (twist() < 0.5) { // comp goes first against its neighbour (order determined randomly)
                    cpay += payoff1st[world[cr*L + cc]*8 + world[cnr*L + cnc]]; // add to comp's payoff
                    if (world[cr*L + cc] < 4) { // comp cooperates first 
                        if ((world[cnr*L + cnc] == 0) || (world[cnr*L + cnc] == 1) || (world[cnr*L + cnc] == 4) || (world[cnr*L + cnc] == 5)) { // comp neighbour responds with cooperation
                            outcome[0]++; // increment outcome 0 (CC)
                        }
                        else { // comp neighbour responds with defection
                            outcome[1]++; // increment outcome 1 (CD)
                        }
                    }
                    else { // comp defects first
                        if ((world[cnr*L + cnc] == 0) || (world[cnr*L + cnc] == 2) || (world[cnr*L + cnc] == 4) || (world[cnr*L + cnc] == 6)) { // comp neighbour responds with cooperation
                            outcome[2]++; // increment outcome 2 (DC)
                        }
                        else { // comp neighbour responds with defection
                            outcome[3]++; // increment outcome 3 (DD)
                        }
                    }
                }
                else { // comp goes second against its neighbour (order determined randomly)
                    cpay += payoff2nd[world[cr*L + cc]*8 + world[cnr*L + cnc]]; // add to comp's payoff
                    if (world[cnr*L + cnc] < 4) { // comp neighbour cooperates first
                        if ((world[cr*L + cc] == 0) || (world[cr*L + cc] == 1) || (world[cr*L + cc] == 4) || (world[cr*L + cc] == 5)) { // comp responds with cooperation
                            outcome[0]++; // increment outcome 0 (CC)
                        }
                        else { // comp responds with defection
                            outcome[1]++; // increment outcome 1 (CD)
                        }    
                    }
                    else { // comp neighbour defects first
                        if ((world[cr*L + cc] == 0) || (world[cr*L + cc] == 2) || (world[cr*L + cc] == 4) || (world[cr*L + cc] == 6)) { // comp responds with cooperation 
                            outcome[2]++; // increment outcome 2 (DC)
                        }
                        else { // comp responds with defection
                            outcome[3]++; // increment outcome 3 (DD)
                        }
                    }
                }
            }
            fpay /= 4.0; // convert fpay from cumulative to average payoff
            cpay /= 4.0; // convert cpay from cumulative to average payoff

            // 2.3.5. replacement
            if (beta < 0.0) { // Replicator rule (set beta to -1 to invoke replicator rule)
                if (alpha > 0.0) { // i.e., if they have all the same payoffs, no replacement occurs
                    if ((cpay > fpay) && (twist() < (cpay - fpay)/alpha)) { // i.e., if comp replaces focal (Replicator rule)
                        abund[world[fr*L + fc]]--; // decrement abundance of focal's strategy
                        abund[world[cr*L + cc]]++; // increment abundance of comp's strategy
                        world[fr*L + fc] = world[cr*L + cc]; // replacement (due to competition) occurs
                    }
                }
            }
            else { // Fermi rule
                if (twist() < 1.0/(1.0 + exp(-beta*(cpay - fpay)))) { // i.e., if comp replaces focal (Fermi rule)
                    abund[world[fr*L + fc]]--; // decrement abundance of focal's strategy
                    abund[world[cr*L + cc]]++; // increment abundance of comp's strategy
                    world[fr*L + fc] = world[cr*L + cc]; // replacement (due to competition) occurs
                }
            }
                        
            // 2.3.6. mutation
            if (twist() < mu) { // i.e., if mutation occurs
                if (seq) { // sequential version of the model (all 8 strategies available)
                    mutant = static_cast<int>(twist()*8); // strategy of mutant (random int between 0 and 7)
                    abund[world[fr*L + fc]]--; // decrement abundance of focal's strategy
                    abund[mutant]++; // increment abundance of mutant's strategy
                    world[fr*L + fc] = mutant; // replacement (due to mutation) occurs
                }
                else { // simultaneous version of the model (only strategies 0 (CC) and 7 (DD) are available)
                    if (twist() < 0.5) { // mutant is strategy 0
                        abund[world[fr*L + fc]]--; // decrement abundance of focal's strategy
                        abund[0]++; // increment abundance of mutant's strategy
                        world[fr*L + fc] = 0; // replacement (due to mutation) occurs
                    }
                    else { // mutant is strategy 7
                        abund[world[fr*L + fc]]--; // decrement abundance of focal's strategy
                        abund[7]++; // increment abundance of mutant's strategy
                        world[fr*L + fc] = 7; // replacement (due to mutation) occurs
                    }
                }
            }
            
        } // end of t loop
        
        // 2.3.6. record time-series data
        for (int j = 0; j < 8; j++) {tsabund[(g + 1)*8 + j] = abund[j];} // add row (g + 1) to time-series-abundance data
        for (int j = 0; j < 4; j++) {tsoutcome[g*4 + j] = outcome[j];} // add row g to time-series-outcome data
        
        // 2.3.7. record snapshots   
        if (snapon && (((g + 1)%snapint) == 0)) {
            snapID++;
            for (int j = 0; j < L*L; j++) {
                snapshot[snapID*L*L + j] = world[j]; // add row snapID (snapshot of current conditions)
            }
        }            
        
    } // end of g loop; end of simulation
        
    // 2.4. OUTPUT RESULTS TO MATLAB
    
    // 2.4.1. output tsabund
    ii = 0;
    for (int j = 0; j < 8; j++) {
        for (int i = 0; i < (gmax + 1); i++) {
            tsabund_lhs[ii] = tsabund[i*8 + j]; // output in column-major order (Matlab format)
            ii++;
        }
    }
    
    // 2.4.2. output tsoutcome
    ii = 0;
    for (int j = 0; j < 4; j++) {
        for (int i = 0; i < gmax; i++) {
            tsoutcome_lhs[ii] = tsoutcome[i*4 + j]; // output in column-major order (Matlab format)
            ii++;
        }
    }
    
    // 2.4.3. output snapshot
    if (snapon) {
        ii = 0;
        for (int j = 0; j < L*L; j++) {
            for (int i = 0; i < nsnap; i++) {
                snapshot_lhs[ii] = snapshot[i*L*L + j]; // output in column-major order (Matlab format)
                ii++;
            }
        }
    }
    else { // dummy output if there are no snapshots
        snapshot_lhs[0] = -1;
    }
        
    
    // 2.5. FREE UP MEMORY THAT IS NO LONGER NEEDED (DELETE POINTERS)
    
    delete[] world;         world = NULL;
    delete[] tsabund;       tsabund = NULL;
    delete[] tsoutcome;     tsoutcome = NULL;
    delete[] snapshot;      snapshot = NULL;

}

// 3. SUB-ROUTINES AND SUPPLEMENTARY FUNCTIONS

// 3.9. twist
double genrand_real2(void); // function prototype
double twist(void) 
// function: returns a random double on [0, 1) using Mersenne Twister
// Mersenne Twister code given in Appendix (section 5)
{
    return genrand_real2();
}

// 4. MEX FUNCTION
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    // arguments
    double *tsabund_lhs, *tsoutcome_lhs, *snapshot_lhs; // lhs arguments
    double *R_rhs, *S_rhs, *T_rhs, *P_rhs, *L_rhs, *startcond_rhs, *beta_rhs, *mu_rhs, *gmax_rhs, *seq_rhs, *sp_rhs, *snapon_rhs, *snapint_rhs; // rhs arguments
    
    // simple error checking
    if (nrhs != 13) {
        mexErrMsgTxt("model requires 13 input arguments");
    }
    else if (nlhs > 3) {
        mexErrMsgTxt("model requires 0-3 output arguments");
    }
    
    // get some useful quantitities
    const int gmax_in = static_cast<int>(*mxGetPr(prhs[8]));
    const int L_in = static_cast<int>(*mxGetPr(prhs[4]));
    const int snapon_in = static_cast<int>(*mxGetPr(prhs[11]));
    const int snapint_in = static_cast<int>(*mxGetPr(prhs[12]));
    const int nsnap_in = static_cast<int>(gmax_in/snapint_in + 1);
    
    // initiate mrows and ncols, used to size output matrices
    int mrows;
    int ncols;
    
    // prepare zeroth output matrix (tsabund_lhs)
    mrows = gmax_in + 1;
    ncols = 8;
    plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL); // tsabund_lhs
    
    // prepare first output matrix (tsoutcome_lhs)
    mrows = gmax_in;
    ncols = 4;
    plhs[1] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL); // tsoutcome_lhs
    
    // prepare second output matrix (snapshot_lhs)
    if (snapon_in) { // record snapshots
        mrows = nsnap_in;
        ncols = L_in*L_in; 
    }
    else { // don't record snapshots
        mrows = 1;
        ncols = 1;
    }
    plhs[2] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL); // snapshot_lhs
    
    // define lhs variables
    tsabund_lhs = mxGetPr(plhs[0]);
    tsoutcome_lhs = mxGetPr(plhs[1]);
    snapshot_lhs = mxGetPr(plhs[2]);
            
    // define rhs variables
    R_rhs = mxGetPr(prhs[0]);
    S_rhs = mxGetPr(prhs[1]);
    T_rhs = mxGetPr(prhs[2]);
    P_rhs = mxGetPr(prhs[3]);
    L_rhs = mxGetPr(prhs[4]);
    startcond_rhs = mxGetPr(prhs[5]);
    beta_rhs = mxGetPr(prhs[6]);
    mu_rhs = mxGetPr(prhs[7]);
    gmax_rhs = mxGetPr(prhs[8]);
    seq_rhs = mxGetPr(prhs[9]);
    sp_rhs = mxGetPr(prhs[10]);
    snapon_rhs = mxGetPr(prhs[11]);
    snapint_rhs = mxGetPr(prhs[12]);    
    
    // run simulation
    sim(tsabund_lhs, tsoutcome_lhs, snapshot_lhs,
            R_rhs, S_rhs, T_rhs, P_rhs, L_rhs, startcond_rhs, beta_rhs, mu_rhs, gmax_rhs, seq_rhs, sp_rhs, snapon_rhs, snapint_rhs); 
}

// 5. APPENDIX - Mersenne Twister code

// Original code by Takuji Nishimura and Makoto Matsumoto (2002)
// Included and modified as allowed in the material below:

/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#include <stdio.h>

// MODIFIED: added to allow for seeding from clock and process ID
#include <ctime>
#include <process.h>
#define GETPID _getpid

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

// OMITTED: "void init_by_array(unsigned long init_key[], int key_length)" 

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(getpid()*time(NULL)); //MODIFIED: Original: init_genrand(5489UL); /* a default initial seed is used */    
            
        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

// OMITTED: "long genrand_int31(void)"
// OMITTED: "genrand_real1(void)"

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

// OMITTED: "genrand_real3(void)"
// OMITTED: "genrand_res53(void)"
// OMITTED: "main(void)"