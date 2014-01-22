#!/bin/bash
# sigma.plt: plot stress tensor components in one or more simulations
# use: sigma.plt mdrun1.r [mdrun2.r ...]
gnuplot -persist <<EOF
plot "<grep -h sigma_xx $*" u 2 w l, "<grep -h sigma_yy $*" u 2 w l, \
"<grep -h sigma_zz $*" u 2 w l, "<grep -h sigma_xy $*" u 2 w l, \
"<grep -h sigma_yz $*" u 2 w l, "<grep -h sigma_xz $*" u 2 w l, 0
fit sxx "<grep -h '<sigma_xx' $*" u 0:2 via sxx
fit syy "<grep -h '<sigma_yy' $*" u 0:2 via syy
fit szz "<grep -h '<sigma_zz' $*" u 0:2 via szz
fit sxy "<grep -h '<sigma_xy' $*" u 0:2 via sxy
fit syz "<grep -h '<sigma_yz' $*" u 0:2 via syz
fit sxz "<grep -h '<sigma_xz' $*" u 0:2 via sxz
print "avg sxx: ", sxx
print "avg syy: ", syy
print "avg szz: ", szz
print "avg sxy: ", sxy
print "avg syz: ", syz
print "avg sxz: ", sxz
EOF
