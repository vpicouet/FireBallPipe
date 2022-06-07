#/bin/bash
dir=~/fireball2-shared/data/dobc/flight/
while inotifywait -e close_write $dir/180922/ ; do
    ~/FireBallPipe/MergeTP.py $dir
done
