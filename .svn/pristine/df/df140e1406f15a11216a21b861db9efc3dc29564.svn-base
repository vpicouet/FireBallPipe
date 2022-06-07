#/bin/bash
#with CIP
file=~/fireball2-shared/data/dobc/flight/1Mbit/rxdata/waterpress_master.csv
# with wifi
#file=~/fireball2-shared/data/dobc/flight/180922/waterpress.csv
while inotifywait -e close_write $file; do
    ~/FireBallPipe/convert_pressure.py $file
done
