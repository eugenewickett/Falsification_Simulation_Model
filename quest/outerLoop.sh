#!/bin/bash
## this next line reads tab-delimited (\t) and expects to get two values per line, 
## saved as P1 and P2
# ./outerLoop.sh [mainInner.sh][replication amount]

for i in seq 1 $2 ;
	do sbatch $1 ;
done
