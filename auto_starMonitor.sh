#/bin/bash

_term() { 
  echo "Caught SIGTERM signal!" 
  pkill -TERM -P "$child1" 2>/dev/null
  pkill -TERM -P "$child2" 2>/dev/null
}

trap _term TERM INT

dir_r12=~/fireball2-shared/data/gobc/ground/hagrid/fbground_r12/profiles/
dir_r13=~/fireball2-shared/data/gobc/ground/hagrid/fbground_r13/profiles/

python ~/FireBallPipe/starMonitor.py $dir_r12 &
child1=$! 
python ~/FireBallPipe/starMonitor.py $dir_r13 &
child2=$!

wait "$child1"
wait "$child2"
#trap - TERM INT
#wait "$child1"
#wait "$child2"

    

