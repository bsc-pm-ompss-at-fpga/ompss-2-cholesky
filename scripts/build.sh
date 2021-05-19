#!/bin/bash -e

if [ "$BOARD" == "" ]; then
  echo "BOARD environment variable not defined"
  exit 1
elif [ "$FPGA_HWRUNTIME" == "" ]; then
  echo "FPGA_HWRUNTIME environment variable not defined"
  exit 1
fi

OUT_DIR=$(pwd -P)/build
RES_FILE=$(pwd -P)/resources_results.json

if [ "$NODE_NAME" == "hca-ikergune" ]; then
  module load openblas/0.3.6_arm64
  #FIXME: Ensure that we are using the arm64 compilers as openblas lib targets this arch
  export CROSS_COMPILE='aarch64-linux-gnu-'
elif [ "$NODE_NAME" == "llebeig" ]; then
  module load openblas/0.3.13/sequential
elif [ "$NODE_NAME" == "xaloc" ]; then
  module load openblas/0.3.10/sequential
else
  echo "NODE_NAME ('$NODE_NAME') not supported. Cannot load OpenBLAS."
  exit 1
fi

mkdir -p $OUT_DIR

if [ "$1" == "binary" ]; then
  #Only build the binaries
  make clean
  make cholesky-p cholesky-d
  mv cholesky-p cholesky-d $OUT_DIR
else
  make clean
  make bitstream-p LDFLAGS=--Wf,--disable_utilization_check
  mv cholesky-p $OUT_DIR
  mv cholesky_ait/cholesky.bin $OUT_DIR/bitstream.bin
  mv cholesky_ait/cholesky.bit $OUT_DIR/bitstream.bit

  printf "{\"version\": \"${FPGA_HWRUNTIME}_${BLOCK_SIZE}\", " >>$RES_FILE
  printf "\"hwruntime\": \"${FPGA_HWRUNTIME}\", " >>$RES_FILE
  printf "\"accels_freq\": \"${FPGA_CLOCK}\", " >>$RES_FILE
  printf "\"memory_port_width\": \"${FPGA_MEMORY_PORT_WIDTH}\", " >>$RES_FILE
  printf "\"benchmark\": \"cholesky\", " >>$RES_FILE
  printf "\"BRAM\": \"" >>$RES_FILE
  grep "BRAM" cholesky_ait/cholesky.resources-hls.txt | awk '{printf $2}' >>$RES_FILE
  printf "\", \"DSP\": \"" >>$RES_FILE
  grep "DSP" cholesky_ait/cholesky.resources-hls.txt | awk '{printf $2}' >>$RES_FILE
  printf "\", \"FF\": \"" >>$RES_FILE
  grep "FF" cholesky_ait/cholesky.resources-hls.txt | awk '{printf $2}' >>$RES_FILE
  printf "\", \"LUT\": \"" >>$RES_FILE
  grep "LUT" cholesky_ait/cholesky.resources-hls.txt | awk '{printf $2}' >>$RES_FILE
  printf "\", \"BRAM_IMPL\": \"" >>$RES_FILE
  grep "BRAM" cholesky_ait/cholesky.resources-impl.txt | awk '{printf $2}' >>$RES_FILE
  printf "\", \"DSP_IMPL\": \"" >>$RES_FILE
  grep "DSP" cholesky_ait/cholesky.resources-impl.txt | awk '{printf $2}' >>$RES_FILE
  printf "\", \"FF_IMPL\": \"" >>$RES_FILE
  grep "FF" cholesky_ait/cholesky.resources-impl.txt | awk '{printf $2}' >>$RES_FILE
  printf "\", \"LUT_IMPL\": \"" >>$RES_FILE
  grep "LUT" cholesky_ait/cholesky.resources-impl.txt | awk '{printf $2}' >>$RES_FILE
  printf "\", \"WNS\": \"" >>$RES_FILE
  grep "WNS" cholesky_ait/cholesky.timing-impl.txt | awk '{printf $2}' >>$RES_FILE
  printf "\", \"TNS\": \"" >>$RES_FILE
  grep "TNS" cholesky_ait/cholesky.timing-impl.txt | awk '{printf $2}' >>$RES_FILE
  printf "\", \"NUM_ENDPOINTS\": \"" >>$RES_FILE
  grep "NUM_ENDPOINTS" cholesky_ait/cholesky.timing-impl.txt | awk '{printf $2}' >>$RES_FILE
  printf "\", \"NUM_FAIL_ENDPOINTS\": \"" >>$RES_FILE
  grep "NUM_FAIL_ENDPOINTS" cholesky_ait/cholesky.timing-impl.txt | awk '{printf $2}' >>$RES_FILE
  printf "\"},\n" >>$RES_FILE
fi
