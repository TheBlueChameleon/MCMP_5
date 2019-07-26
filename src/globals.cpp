/* globals.hpp
 * behavioural parameters and global macros, debug file and GSL RNG
 */

// ========================================================================= //
// dependencies

#include <iostream>
#include <ctime>
#include <string>
#include <fstream>
#include <cmath>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "globals.hpp"

// ========================================================================= //
// global variables

// ------------------------------------------------------------------------- //
// nature constants

const double PI = std::atan(1.0) * 4;
const std::complex<double> I = std::complex<double>(0, 1);

// ------------------------------------------------------------------------- //
// simulation behaviour

extern const double outerT_lo = 0.30;
extern const double outerT_hi = 4.00;
extern const double outerT_dT = 0.10;

extern const double innerT_lo = 0.60;
extern const double innerT_hi = 1.40;
extern const double innerT_dT = 0.02;

// ------------------------------------------------------------------------- //
// log behaviour

const std::string   REPORT_DIR     = "./txtout/";
const std::string   DEBUG_FILENAME = "MCMP_4_Debug.txt";

std::ofstream hDebug(DEBUG_FILENAME);

// ------------------------------------------------------------------------- //
// GSL RNG

bool      RNG_initialized = false;
gsl_rng * RNG;


// ========================================================================= //
// proc

void RNG_quit ();

// ------------------------------------------------------------------------- //

void RNG_init () {
  gsl_rng_env_setup();
  RNG = gsl_rng_alloc (gsl_rng_default);    // default is MT19937
  
  time_t seed = std::time(nullptr);
  gsl_rng_set(RNG, seed);
  hDebug << "RNG of type " << gsl_rng_name (RNG) << " initialized with seed " << seed << std::endl;
  
  std::atexit(RNG_quit);
  
  RNG_initialized = true;
}

// ------------------------------------------------------------------------- //

void RNG_quit () {
  gsl_rng_free (RNG);
  RNG_initialized = false;
}
