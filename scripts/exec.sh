#!/bin/bash -el

RES_FILE=$(pwd -P)/test_results.json

for EXEC_MODE in d p; do
  for MATRIX_SIZE in 5120 10240; do
    echo "=== Check mode: ${EXEC_MODE}, msize: ${MATRIX_SIZE} ==="
    ##NOTE: Check == 2 -> enables the warm-up mode
    CHECK=$([ "$MATRIX_SIZE" == "5120" ] && echo 1 || echo 2)
	NX_ARGS="--summary --smp-workers=12 --fpga-alloc-pool-size=4G" timeout --preserve-status 150s ./build/cholesky-${EXEC_MODE} ${MATRIX_SIZE} ${CHECK}
    cat test_result.json >>$RES_FILE
    echo "," >>$RES_FILE
  done
done
