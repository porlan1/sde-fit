/* sdef_NLL.h */
#ifndef SDEF_NLL_H
#define SDEF_NLL_H

#include "sdef_misc.h"
#include "sdef_Fokker_Planck.h"

double NLLF_FP(Model_data MD, parameters par);
double NLLF_MBB(Model_data MD);

#endif /* sdef_NLL.h */