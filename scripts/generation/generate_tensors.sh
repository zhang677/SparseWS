FOLDER="/home/zgh23/code/SparseWS/data/origin/tensors/"
# K=38955429
# L=38955429
# I=532
# k=7791085 # 7791085 # 3895543 # 
# l=7791085 # 7791085 # 3895543 #
K=1504
L=42390
k=300
l=8479
J=$1 # 4 # 4 # 8
j=$2 # 2 # 1 # 3
Tname=facebook
Dname="/scratch/zgh23/sparse_ten/D-$Tname-$J.txt"
Cname="/scratch/zgh23/sparse_ten/C-$Tname-$J.txt"
Dmtx="/scratch/zgh23/sparse_ten/D-$Tname-$J.mtx"
Cmtx="/scratch/zgh23/sparse_ten/C-$Tname-$J.mtx"
./generate_rec $K $J $k $j $FOLDER 0
cp ${FOLDER}rec\_$k\_$K\_$j\_$J\_0.mtx $Cname
# Remove the first 2 lines 
sed -i '1,2d' $Cname
cp ${FOLDER}rec\_$k\_$K\_$j\_$J\_0.mtx $Cmtx
./generate_rec $L $J $l $j $FOLDER 42
cp ${FOLDER}rec\_$l\_$L\_$j\_$J\_42.mtx $Dname
sed -i '1,2d' $Dname
cp ${FOLDER}rec\_$l\_$L\_$j\_$J\_42.mtx $Dmtx
# ./generate_rec $I $J $i $j $FOLDER 42
# cp ${FOLDER}rec\_$i\_$I\_$j\_$J\_42.mtx $Aname
# sed -i '1,2d' $Aname
# Copy the generated tensors to the scratch directory


