#!/usr/bin/env bash

usage(){
cat << EOF 

     usage:     bash mod_init.sh -i /path/to/flu/init.sh 

EOF
}

while [ "$1" != "" ]; do
    case $1 in
        -i | --input )
            shift
            input=$1 
            ;;           
        -h | --help )    usage
            exit
        ;;
        * )              usage
            exit 1
    esac
    shift
done

echo $input

#backup init.sh
cp $input $input.backup

# settings to change
SINGLE_LOCAL_PROC=10    # local maximum processes
MIN_LEN=75          # minimum read length
MIN_FI=0.05     # minimum insertion variant frequency
MIN_FD=0.05     # minimum deletion variant frequency
MIN_F=0.08      # minimum frequency for single nucleotide variants
let "DOUBLE_LOCAL_PROC = $SINGLE_LOCAL_PROC /2"

# find vars
SLP=$(grep -irn "SINGLE_LOCAL_PROC=" $input | tr ":" "\t" | cut -f 1)
SLP2=$(grep -irn "DOUBLE_LOCAL_PROC=" $input | tr ":" "\t" | cut -f 1)
LEN=$(grep -irn "MIN_LEN=" $input | tr ":" "\t" | cut -f 1)
MFI=$(grep -irn "MIN_FI=" $input | tr ":" "\t" | cut -f 1)
MFD=$(grep -irn "MIN_FD=" $input | tr ":" "\t" | cut -f 1)
MF=$(grep -irn "MIN_F=" $input | tr ":" "\t" | cut -f 1)

# replace with these settings
sed -i "$SLP s/.*/SINGLE_LOCAL_PROC=$SINGLE_LOCAL_PROC/" $input
sed -i "$SLP2 s/.*/DOUBLE_LOCAL_PROC=$DOUBLE_LOCAL_PROC/" $input
sed -i "$LEN s/.*/MIN_LEN=$MIN_LEN/" $input
sed -i "$MFI s/.*/MIN_FI=$MIN_FI/" $input
sed -i "$MFD s/.*/MIN_FD=$MIN_FD/" $input
sed -i "$MF s/.*/MIN_F=$MIN_F/" $input

# insert these new ones
sed -i "$ a MIN_CONS_SUPPORT=20" $input
sed -i "$ a MIN_CONS_QUALITY=15" $input
sed -i '$ a ALIGN_AMENDED=1' $input

