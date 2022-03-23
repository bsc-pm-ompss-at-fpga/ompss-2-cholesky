#!/bin/bash -e

if [ "$BOARD" == "" ]; then
  echo "BOARD environment variable not defined"
  exit 1
elif [ "$FPGA_HWRUNTIME" == "" ]; then
  echo "FPGA_HWRUNTIME environment variable not defined"
  exit 1
fi

PROG_NAME=cholesky
OUT_DIR=$(pwd -P)/build
RES_FILE=$(pwd -P)/resources_results.json

if [ "$NODE_NAME" == "hca-ikergune" ]; then
  module load openblas/0.3.6_arm64
  #FIXME: Ensure that we are using the arm64 compilers as openblas lib targets this arch
  export CROSS_COMPILE='aarch64-linux-gnu-'
elif [ "$NODE_NAME" == "llebeig" ]; then
  module load openblas/0.3.13/sequential
elif [ "$NODE_NAME" == "xaloc" ]; then
  module load openblas/0.3.12/sequential
elif [ "$NODE_NAME" == "quar" ]; then
  module load openblas/0.3.17/sequential
else
  echo "NODE_NAME ('$NODE_NAME') not supported. Cannot load OpenBLAS."
  exit 1
fi

mkdir -p $OUT_DIR

if [ "$1" == "binary" ]; then
  #Only build the binaries
  make clean
  make ${PROG_NAME}-p ${PROG_NAME}-d
  mv ${PROG_NAME}-p ${PROG_NAME}-d $OUT_DIR
elif [ "$1" == "design" ]; then
  #Only generate the design
  make clean
  make design-p design-i

  #Remove OUT_DIR directory since we are not generating output products
  rm -rf $OUT_DIR
else
  make clean
  make bitstream-p LDFLAGS=--Wf,--disable_utilization_check

  mv ${PROG_NAME}_ait/${PROG_NAME}.bin $OUT_DIR/bitstream.bin
  mv ${PROG_NAME}_ait/${PROG_NAME}.bit $OUT_DIR/bitstream.bit
  mv ${PROG_NAME}_ait/${PROG_NAME}.xtasks.config $OUT_DIR/xtasks.config

  printf "{\"version\": \"${FPGA_HWRUNTIME}_${BLOCK_SIZE}\", " >>$RES_FILE
  printf "\"hwruntime\": \"${FPGA_HWRUNTIME}\", " >>$RES_FILE
  printf "\"accels_freq\": \"${FPGA_CLOCK}\", " >>$RES_FILE
  printf "\"memory_port_width\": \"${FPGA_MEMORY_PORT_WIDTH}\", " >>$RES_FILE
  printf "\"benchmark\": \"${PROG_NAME}" >>$RES_FILE
  for PARAM in BRAM DSP FF LUT; do
    printf "\", \"${PARAM}\": \"" >>$RES_FILE
    grep "$PARAM" ${PROG_NAME}_ait/${PROG_NAME}.resources-hls.txt | awk '{printf $2}' >>$RES_FILE
  done
  for PARAM in BRAM DSP FF LUT; do
    printf "\", \"${PARAM}_IMPL\": \"" >>$RES_FILE
    grep "${PARAM}" ${PROG_NAME}_ait/${PROG_NAME}.resources-impl.txt | awk '{printf $2}' >>$RES_FILE
  done
  for PARAM in WNS TNS NUM_ENDPOINTS NUM_FAIL_ENDPOINTS; do
    printf "\", \"${PARAM}\": \"" >>$RES_FILE
    grep "${PARAM}" ${PROG_NAME}_ait/${PROG_NAME}.timing-impl.txt | awk '{printf $2}' >>$RES_FILE
  done
  printf "\"},\n" >>$RES_FILE
fi
