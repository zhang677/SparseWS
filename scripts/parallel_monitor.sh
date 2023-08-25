rt=/home/nfs_data/zhanggh/SparseWS/data/monitor
pids=$(ps -ef | grep test-hash | awk '{print $2}')
readarray -t y <<<"$pids"
arraylength=${#y[@]}-1

# Use for loop to read all pids, except the last (grep process)
# '&' means non-blocking mode, so all processes will run concurrently
# '>>' means append to the file. Record all pids to pids.txt
# The minimal interval is 0.02s
for (( i=0; i<${arraylength}; i++ ));
do
  echo "${y[$i]}" >> $rt/pids.txt 
  psrecord --plot $rt/${y[$i]}.png ${y[$i]} &
done