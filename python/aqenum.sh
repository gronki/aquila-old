#!/bin/bash

set -e
i=1

if [[ $# -ne 2 ]]; then
    echo "usage: cat list.txt | aqenum [basename] [extension]"
    exit
fi

while read fn
do
    ln -srfv "$fn" "$1$i.$2"
    i=$(( i + 1 ))
done
