#pragma once
#include "Rtypes.h"

/**
 * Contains daq data structures
 */

#define CAEN_1742_GR 4
#define CAEN_1742_CH 32
#define CAEN_1742_LN 1024

struct caen_1742 {
  ULong64_t system_clock;
  ULong64_t device_clock[CAEN_1742_CH];
  UShort_t trace[CAEN_1742_CH][CAEN_1742_LN];
  UShort_t trigger[CAEN_1742_GR][CAEN_1742_LN];
};