rt=/home/nfs_data/zhanggh/SparseWS/data/monitor
pids=$(ps -ef | grep test-hash | awk '{print $2}')
readarray -t y <<<"$pids"
arraylength=${#y[@]}-1

# use for loop to read all values and indexes
for (( i=0; i<${arraylength}; i++ ));
do
  echo "${y[$i]}" >> $rt/pids.txt 
  psrecord --plot $rt/${y[$i]}.png --log $rt/${y[$i]}.txt ${y[$i]} &
done