#!/bin/bash -ex
module load xl
INSTALL_DIR=/gpfs/lb/provisioned/FEPS/shared/installXL/
if [ -e feaprog.o ]; then
  rm feaprog.o
fi
mpixlcxx -DHAVE_CONFIG_H -D_HAVE_ZOLTAN_ -D_HAVE_PARMETIS_ -DFMDB_PARALLEL -DMESHMODEL -DFMDB \
         -I${INSTALL_DIR}/include \
         -g -c feaprog.cc -o feaprog.o 
mpixlcxx  feaprog.o -L/bgsys/drivers/ppcfloor/comm/xl/lib -L${INSTALL_DIR}/lib/ \
         -lFMDB -lSCORECModel -lmeshModel -lSCORECUtil -lzoltan -lipcomman -lparmetis -lmetis \
         -o feaprog
