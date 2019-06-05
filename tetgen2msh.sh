#!/bin/sh

if [ "X$1" == "X" ]; then
	echo "no mesh specified"
	return 1;
fi

# Strip common suffixes
BASE=${1%%.node}
BASE=${BASE%%.ele}
BASE=${BASE%%.msh}

ELE=$BASE.ele
NODE=$BASE.node
MSH=$BASE.msh

if [ ! -f $ELE ]; then
	echo "$ELE does not exist"
	return 1;
fi

if [ ! -f $NODE ]; then
	echo "$NODE does not exist"
	return 1;
fi

# Parse conductivities
C1="0.33"
C2="0.0042"
C3="1.79"
C4="0.33"

while [ ! -z "$2" ]; do
	v=${2:2}
	case $2 in
		(1=*) C1=$v;;
		(2=*) C2=$v;;
		(3=*) C3=$v;;
		(4=*) C4=$v;;
	esac
	shift
done

echo "Conductivities: 1=${C1}, 2=${C2}, 3=${C3}, 4=${C4}"
echo Generating nodes

sed 's/#.*//' $NODE | awk \
	'BEGIN { nn = 0 }
	nn != 0 && NF == 4 {
		print $1,$2,$3,$4
		next
	}
	nn == 0 && NF == 4 && $2 == 3 {
		nn = $1
		print nn
	}' > $MSH

echo Generating elements

sed 's/#.*//' $ELE | awk \
	'BEGIN { nn = 0 }
	nn != 0 && NF == 6 {
		print $1,$3,$2,$4,$5
		next
	}
	nn == 0 && NF == 3 && $2 == 4 {
		nn = $1
		print nn,4
	}' >> $MSH

echo Generating conductivity

sed 's/#.*//' $ELE | awk \
	'BEGIN { nn = 0 }
	nn != 0 && NF == 6 {
		if (int($6) == 1) print $1,"'$C1'"
		if (int($6) == 2) print $1,"'$C2'"
		if (int($6) == 3) print $1,"'$C3'"
		if (int($6) == 4) print $1,"'$C4'"
		next
	}
	nn == 0 && NF == 3 && $2 == 4 {
		nn = $1
		print nn
	}' >> $MSH
