WORKDIR=`pwd`
EXECP=${CAMPARI_HOME}/bin/${ARCH}

mkdir ${WORKDIR}/RESULTS

for i in 20 30 50 75 100 150 250;
do
  dn=`printf %03d $i`
  cd ${WORKDIR}
  rm ${WORKDIR}/leu.in 
  echo "ACE" > ${WORKDIR}/leu.in
  for j in `seq 1 $i`;
  do
    echo "LEU" >> ${WORKDIR}/leu.in
  done
  echo -e "NME\nEND" >> ${WORKDIR}/leu.in
  ${EXECP}/campari -k ${WORKDIR}/IPPChain.key >> ${WORKDIR}/IPPChain.log
  for j in KRATKY PERSISTENCE POLYAVG INTSCAL RETEHIST
  do
    mv ${j}.dat ${WORKDIR}/RESULTS/${j}_${dn}.dat
  done
done

cd ${WORKDIR}
rm *.dat *.tmp *.int *.pdb

for i in 20 30 50 75 100 150 250;
do
  dn=`printf %03d $i`
  val=`grep -v -F "#" ${WORKDIR}/RESULTS/POLYAVG_${dn}.dat | awk '{print $1}'`
  echo "$i  $val" >> IPPscaling.dat
done
