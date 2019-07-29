#ifndef __CHOLESKU_FPGA_H__
#define __CHOLESKU_FPGA_H__

const int ts = 64; // tile size

#if defined(USE_DOUBLE)
#  define type_t     double
#  define ELEM_T_STR "double"
#else
#  define type_t     float
#  define ELEM_T_STR "float"
#endif

#endif /* __CHOLESKU_FPGA_H__ */
