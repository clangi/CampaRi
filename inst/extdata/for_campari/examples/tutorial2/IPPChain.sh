WORKDIR=`pwd`
CAMPARIEXE=${CAMPARI_HOME}/bin/${ARCH}/campari_threads
AMINO=LEU

mkdir ${WORKDIR}/RESULTS

for i in 20 60 120 200 300 420 560 720 900;
do
  dn=`printf %03d $i`
  cd ${WORKDIR}
  rm ${WORKDIR}/tutorial2.in 
  echo "${AMINO}_N" > ${WORKDIR}/tutorial2.in
  for j in `seq 3 $i`; # subtract two terminal residues from loop
  do
    echo "${AMINO}" >> ${WORKDIR}/tutorial2.in
  done
  echo -e "${AMINO}_C\nEND" >> ${WORKDIR}/tutorial2.in
  ${CAMPARIEXE} -k ${WORKDIR}/tutorial2.key >> ${WORKDIR}/tutorial2.log
  for j in DENSPROF PERSISTENCE POLYAVG INTSCAL RETEHIST
  do
    mv ${j}.dat ${WORKDIR}/RESULTS/${j}_${dn}.dat
  done
done

cd ${WORKDIR}
rm *.dat *.tmp *.int *.pdb

for i in 20 60 120 200 300 420 560 720 900;
do
  dn=`printf %03d $i`
  val=`grep -v -F "#" ${WORKDIR}/RESULTS/POLYAVG_${dn}.dat | awk '{print $1}'`
  echo "$i  $val" >> IPPscaling.dat
done
