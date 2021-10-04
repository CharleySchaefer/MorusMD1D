/*
  pname: name of the program, argv[0] to main
*/
#include "printHelp.h"
void printHelp(char *pname) {
    printf("\n  USAGE: %s <arguments>\n", pname);
    printf("\n  ARGUMENTS:\n");
    printf("  --help\n");
    printf("  --version\n");
    printf("\n  INITIALISATION\n");
    printf("  --init-continue (to be implemented)\n");
    printf("  --rseed <int>: initialise random number generator with defined seed.\n");  
    printf("\n  EXPORTING\n");
    printf("  simulation data is exported to \'timeprogress.out\' at a given path.\n");
    printf("  --outdir <path>: default=\'.\'");
    printf("  --Ntime <int>:   number of time steps\n");
    printf("  --Nprint <int>:  export to timeprogress after this number of time steps.\n");
    printf("  --export-logt (overrides --Nprint): time interval exporting evenly speced on a log scale\n");
    printf("  --export-sequence: export coordinate of each bead in timeprogress.out\n");
    printf("\n  CALCULATING mu(t): mu(t)determines relaxation rate in linear viscoelasticity\n");
    printf("    mu(t) will be calculated on a time interval [tL, tU], with Nbins bins.\n");
    printf("    the bin values are updated after Nsteps simulation steps.\n");
    printf("  --export-mu:     turn on mu(t) tracking\n");
    printf("  --mu-time-steps <Nbins (integer)>.\n");
    printf("  --mu-time-L <tL(float)>.\n");
    printf("  --mu-time-U <tU(float)>\n");
    printf("  --mu-Nbins  <Nbins(integer)>\n");
    printf("\n  CHAIN PROPERTIES\n");
    printf("  A chain of N beads represents a chain of Ze entanglements with Zs stickers\n");
    printf("  in an arbitrary sequence (of length N). The exact representation may be varied\n");
    printf("  by the boundary condition (documentation in the code).\n");
    printf("  --Ze <int>            ; number of entanglements per chain \n");
    printf("  --sequence <string with zeros and ones>: 0=non-sticky bead, 1=sticker\n");
    printf("  --lambda-max  <value> ; maximum chain extension ratio (default infinity)\n");
    printf("  --set-boundary <char>: boundary A, B(default) or C.\n");
    printf("                         see data/demo_extensional_flow.sh\n");
    printf("\n  STICKER PROPERTIES\n");
    printf("  --sticker-sync:        stickers synchronised: either all open or all closed.\n");
    printf("  --sticker-p <float>: fraction of closed stickers.\n");
    printf("  --sticker-kdiss <float>: sticker dissociation time (units of taue).\n");
    printf("  --sticker-kbondswap <float>: rate of sticker bondswapping (units of taue).\n");
    printf("\n  PHYSICAL CONDITIONS\n");
    printf("  --exclude-fluctuations\n");
    printf("  --Wi <float>\n");
    printf("    Weissenberg number for extensional flow \n");
}
