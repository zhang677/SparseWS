#!/bin/bash

cd $1

wget https://s3.us-east-2.amazonaws.com/frostt/frostt_data/nell/nell-1.tns.gz
gunzip nell-1.tns.gz
wget https://s3.us-east-2.amazonaws.com/frostt/frostt_data/flickr/flickr-3d.tns.gz
gunzip flickr-3d.tns.gz
wget https://s3.us-east-2.amazonaws.com/frostt/frostt_data/delicious/delicious-3d.tns.gz
gunzip delicious-3d.tns.gz
wget https://s3.us-east-2.amazonaws.com/frostt/frostt_data/patents/patents.tns.gz
gunzip patents.tns.gz