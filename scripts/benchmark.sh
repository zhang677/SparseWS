METHOD=$1

case $METHOD in
  hash)

  for cap in 262144 524288 1048576; do
    echo $cap
    ../test-hash ../data/origin/2D_27628_bjtcai.mtx ../data/shift/2D_27628_bjtcai-st.mtx 10 $cap 0
  done
  ;;

  coord)
  for cap in 262144 524288 1048576; do
    echo $cap
    ../test-coord ../data/origin/2D_27628_bjtcai.mtx ../data/shift/2D_27628_bjtcai-st.mtx 10 $cap 0
  done
  ;;

esac  