#! /bin/sh
#
# this script generates PDB format files
# for the moving molecule visualization
#
# for each input line generates the PDB
# file from the moving molecule PDB file.
#
# input line fields:
#       $2, $3, $4 - Cartesian coordinates;
#       $5, $6, $7 - Eulerian angles
#
# the input line number is appended to the
# output file name
#
# the moving molecule file must be in
# the standard PDB format
#
# arguments:
#               <1> - the moving molecule file 
#               <2> - the base name for output PDB files
#               <3> - the file with energies and rigid body coordinates
#
# usage: pdbgen <1> <2> <3>
#

OS = `uname`
echo $OS

if [ $# -ne 4 ] ; then
        echo ""
        echo "$0 : Incorrect number of arguments."
        echo ""
        echo " usage: runfade moving.pdb still.pdb baseout configs.list"
#	exit 1
fi
#
echo ""
echo "Generating PDB files in background mode ..."
echo ""
#
awk -v mm=$1 -v sm=$2 -v out=$3 '{ 
       s = $2 " " $3 " " $4 " " $5 " " $6 " " $7 ;
	z="";
        if(NR<1000) z = "0";
	if (NR<100) z = "00";
	if (NR<10) z = "000";
	system( "../exec/pdbshift.$OS " mm " " out z NR ".pdb " s) 
	system( "../exec/FADE.$OS " sm " " out z NR ".pdb " )
}' $4 

