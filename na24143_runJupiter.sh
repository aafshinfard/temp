#!/bin/sh
export PATH="/projects/btl/lcoombe/miniconda3/envs/jupiter/bin:$PATH"

./jupiter ng=95 maxScaff=40 t=16 name=ref.VS.na24143.ONT.shasta.ng95.maxScaff40 ref=grch38.fa fa=na24143.ONT.shasta.fa
./jupiter ng=95 maxScaff=40 t=16 name=ref.VS.na24143.ONT.shasta.arcs.ng95.maxScaff40 ref=grch38.fa fa=na24143.ONT.shasta.arcs.fa
./jupiter ng=95 maxScaff=40 t=16 name=ref.VS.na24143.ONT.shasta.physlr.ng95.maxScaff40 ref=grch38.fa fa=na24143.ONT.shasta.physlr.fa

./jupiter ng=95 maxScaff=40 t=16 name=ref.VS.na24143.pacbio.ng95.maxScaff40 ref=grch38.fa fa=na24143.pacbio.fa
./jupiter ng=95 maxScaff=40 t=16 name=ref.VS.na24143.pacbio.arcs.ng95.maxScaff40 ref=grch38.fa fa=na24143.pacbio.arcs.fa
./jupiter ng=95 maxScaff=40 t=16 name=ref.VS.na24143.pacbio.physlr.ng95.maxScaff40 ref=grch38.fa fa=na24143.pacbio.physlr.fa

./jupiter ng=95 maxScaff=40 t=16 name=ref.VS.na24143.pe_mpet.abyss.ng95.maxScaff40 ref=grch38.fa fa=na24143.pe_mpet.abyss.fa
./jupiter ng=95 maxScaff=40 t=16 name=ref.VS.na24143.pe_mpet.abyss.arcs.ng95.maxScaff40 ref=grch38.fa fa=na24143.pe_mpet.abyss.arcs.fa
./jupiter ng=95 maxScaff=40 t=16 name=ref.VS.na24143.pe_mpet.abyss.physlr.ng95.maxScaff40 ref=grch38.fa fa=na24143.pe_mpet.abyss.physlr.fa

./jupiter ng=95 maxScaff=40 t=16 name=ref.VS.na24143.stLFR.supernova.ng95.maxScaff40 ref=grch38.fa fa=na24143.stLFR.supernova.fa
./jupiter ng=95 maxScaff=40 t=16 name=ref.VS.na24143.stLFR.supernova.physlr.ng95.maxScaff40 ref=grch38.fa fa=na24143.stLFR.supernova.physlr.fa
