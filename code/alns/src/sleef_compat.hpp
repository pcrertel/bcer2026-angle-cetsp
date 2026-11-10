// sleef_compat.hpp
// Use ONLY when your SLEEF does not export the plain scalar alias
// Sleef_atan2d1_u10 (e.g., SLEEF 3.5.x of Ubuntu, which only has the kernels
// with a suffix). Redirects the name to the available scalar kernel, WITHOUT
// changing max_angle_chain.hpp. Pass it with:  -include sleef_compat.hpp
#include <sleef.h>
#if defined(MAXANGLE_USE_SLEEF)
#  if defined(__FMA__)
#    define Sleef_atan2d1_u10 Sleef_atan2d1_u10purecfma
#  else
#    define Sleef_atan2d1_u10 Sleef_atan2d1_u10purec
#  endif
#endif
