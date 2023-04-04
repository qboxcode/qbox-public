#!/bin/bash
# extract <dipole_total> from a run and correct for discontinuities
# use: get_dipole.sh file1.r .. fileN.r

grep -h dipole_total ${*} |awk '{print $2,$3,$4}' > d$$.tmp

alat=$(grep -m 1 '<unit_cell_a_norm>' $1 | awk '{print $2}')
blat=$(grep -m 1 '<unit_cell_b_norm>' $1 | awk '{print $2}')
clat=$(grep -m 1 '<unit_cell_c_norm>' $1 | awk '{print $2}')

# echo $alat $blat $clat

cut -f 1 -d ' ' d$$.tmp | \
  awk -v a=$alat 'NR==1 {printf("%13.10f\n",$1); dm=$1; offset=0}
  NR>1 {d=$1;
      if (d-dm >  a/2) { offset=offset-a; };
      if (d-dm < -a/2) { offset=offset+a; };
      d = d + offset;
      printf("%13.10f\n",d);
      dm=$1}' - > dx$$.tmp

cut -f 2 -d ' ' d$$.tmp | \
  awk -v a=$blat 'NR==1 {printf("%13.10f\n",$1); dm=$1; offset=0}
  NR>1 {d=$1;
      if (d-dm >  a/2) { offset=offset-a; };
      if (d-dm < -a/2) { offset=offset+a; };
      d = d + offset;
      printf("%13.10f\n",d);
      dm=$1}' - > dy$$.tmp

cut -f 3 -d ' ' d$$.tmp | \
  awk -v a=$clat 'NR==1 {printf("%13.10f\n",$1); dm=$1; offset=0}
  NR>1 {d=$1;
      if (d-dm >  a/2) { offset=offset-a; };
      if (d-dm < -a/2) { offset=offset+a; };
      d = d + offset;
      printf("%13.10f\n",d);
      dm=$1}' - > dz$$.tmp

paste dx$$.tmp dy$$.tmp dz$$.tmp
rm d$$.tmp d[xyz]$$.tmp
