#ifndef DAQ_FAST_CORE_INCLUDE_COMMON_HH_
#define DAQ_FAST_CORE_INCLUDE_COMMON_HH_

/*===========================================================================*\

author: Matthias W. Smith
email:  mwsmith2@uw.edu
file:   common.hh

about:  Contains the data structures for several hardware devices in a single
        location.  The header should be included in any program that aims
        to interface with (read or write) data with this daq.

\*===========================================================================*/

#define SIS_3350_CH 4
#define SIS_3350_LN 1024

#define SIS_3302_CH 8
#define SIS_3302_LN 100000

#define SIS_3316_CH 16
#define SIS_3316_GR 4
#define SIS_3316_LN 100000

#define CAEN_1785_CH 8

#define CAEN_6742_GR 2
#define CAEN_6742_CH 18
#define CAEN_6742_LN 1024

#define CAEN_1742_GR 4
#define CAEN_1742_CH 32
#define CAEN_1742_LN 1024

#define CAEN_5720_CH 4
#define CAEN_5720_LN 1024

#define CAEN_5730_CH 8
#define CAEN_5730_LN 500

#define DRS4_CH 4
#define DRS4_LN 1024

#define TEK_SCOPE_CH 4
#define TEK_SCOPE_LN 10000

// NMR specific stuff
#define NMR_FID_LN SIS_3302_LN
#define SHORT_FID_LN 10000
#define SHIM_PLATFORM_CH 28
#define SHIM_FIXED_CH 4
#define RUN_TROLLEY_CH 17
#define RUN_FIXED_CH 378

//--- std includes ----------------------------------------------------------//
#include <vector>
#include <array>
#include <mutex>
#include <cstdarg>
#include <sys/time.h>

//--- other includes --------------------------------------------------------//
#include <boost/variant.hpp>
#include <zmq.hpp>
#include "TFile.h"

//--- projects includes -----------------------------------------------------//
// #include "worker_base.hh"

namespace daq {

// Basic structs
struct test_struct {
  ULong64_t system_clock;
  Double_t  value;
};

struct sis_3350 {
  ULong64_t system_clock;
  ULong64_t device_clock[SIS_3350_CH];
  UShort_t trace[SIS_3350_CH][SIS_3350_LN];
};

struct sis_3302 {
  ULong64_t system_clock;
  ULong64_t device_clock[SIS_3302_CH];
  UShort_t trace[SIS_3302_CH][SIS_3302_LN];
};

struct sis_3316 {
  ULong64_t system_clock;
  ULong64_t device_clock[SIS_3316_CH];
  UShort_t trace[SIS_3316_CH][SIS_3316_LN];
};

struct caen_1785 {
  ULong64_t system_clock;
  ULong64_t device_clock[CAEN_1785_CH];
  UShort_t value[CAEN_1785_CH];
};

struct caen_6742 {
  ULong64_t system_clock;
  ULong64_t device_clock[CAEN_6742_CH];
  UShort_t trace[CAEN_6742_CH][CAEN_6742_LN];
};

struct caen_1742 {
  ULong64_t system_clock;
  ULong64_t device_clock[CAEN_1742_CH];
  UShort_t trace[CAEN_1742_CH][CAEN_1742_LN];
  UShort_t trigger[CAEN_1742_GR][CAEN_1742_LN];
};

struct caen_5720 {
  ULong64_t event_index;
  ULong64_t system_clock;
  UShort_t trace[CAEN_5720_CH][CAEN_5720_LN];
};

struct caen_5730 {
  ULong64_t event_index;
  ULong64_t system_clock;
  UShort_t trace[CAEN_5730_CH][CAEN_5730_LN];
};

struct drs4 {
  ULong64_t system_clock;
  ULong64_t device_clock[DRS4_CH];
  UShort_t trace[DRS4_CH][DRS4_LN];
};

// Built from basic structs
struct event_data {
  std::vector<sis_3350> sis_3350_vec;
  std::vector<sis_3302> sis_3302_vec;
  std::vector<caen_1785> caen_1785_vec;
  std::vector<caen_6742> caen_6742_vec;
  std::vector<caen_1742> caen_1742_vec;
  std::vector<drs4> drs4_vec;
  std::vector<sis_3316> sis_3316_vec;
  std::vector<caen_5720> caen_5720_vec;
  std::vector<caen_5730> caen_5730_vec;
};

// Typedef for all workers - needed by in WorkerList
// typedef boost::variant<WorkerBase<sis_3350> *,
//                        WorkerBase<sis_3302> *,
//                        WorkerBase<caen_1785> *,
//                        WorkerBase<caen_6742> *,
//                        WorkerBase<drs4> *,
//                        WorkerBase<caen_1742> *,
//                        WorkerBase<sis_3316> *,
// 		       WorkerBase<caen_5720> *,
// 		       WorkerBase<caen_5730> *>
// worker_ptr_types;

// A useful define guard for I/O with the vme bus.
extern int vme_dev;
extern std::string vme_path;
extern std::mutex vme_mutex;

// Create a variable for a config directory.
extern std::string conf_dir;

// Set sleep times for data polling threads.
const int short_sleep = 10;
const int long_sleep = 100;
const double sample_period = 0.0001; // in milliseconds

// Set up a global zmq context
extern zmq::context_t msg_context;

inline void light_sleep() {
 usleep(200); // in usec
}

inline void heavy_sleep() {
  usleep(10000); // in usec
}

} // ::daq

#endif
