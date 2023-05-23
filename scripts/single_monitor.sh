rt=/home/nfs_data/zhanggh/SparseWS/data/monitor
pids=$(ps -ef | grep test-hash | awk '{print $2}')
readarray -t y <<<"$pids"
pid=${y[0]}
echo $pid
psrecord --plot $rt/$pid.png --log $rt/$pid.txt $pid