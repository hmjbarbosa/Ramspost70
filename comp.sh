# ***********************************************************************
# Install script Ramspost 5.0
#
# Before run read 'README' file
#
# Use one this options to compile:
#
# - Intel   : ./comp.sh intel
# - G95     : ./comp.sh g95
# - Portand : ./comp.sh pgi
# - Absoft  : ./comp.sh af95
#
# Cleanning generated files (*.o *.mod *.a and binaries):
#
# - Intel   :   ./comp.sh intel clean
# - g95     :   ./comp.sh g95 clean
# - Portand :   ./comp.sh pgi clean
# - Absoft  :   ./comp.sh af95 clean
#
# by Daniel Merli Lamosa PAD/CPTEC/INPE
# 18/05/2005
# **********************************************************************

function Display_not-arg() {

echo "Set a compiler!"
echo ""
echo "Compiling:"
echo " - Intel      : ./comp.sh intel"
echo " - G95        : ./comp.sh g95"
echo " - Portland   : ./comp.sh pgi"
echo " - Absoft     : ./comp.sh af95"
echo ""
echo "Cleannig:"
echo " - Intel      : ./comp.sh intel clean"
echo " - G95        : ./comp.sh g95 clean"
echo " - Portland   : ./comp.sh pgi clean"
echo " - Absoft     : ./comp.sh af95 clean"
echo ""

}


if [ "$1" = "" ]; then
  Display_not-arg
else
  # Argumentos
  ARG1=$1
  ARG2=$2

  cd ./LIB
  make -f Make.utils.opt.${ARG1} $ARG2
  cd ..
  make -f Makefile_60.${ARG1} $ARG2
fi
