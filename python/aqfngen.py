#!/bin/python

from os import read, path, write

for fn in read(0,10000).split('\n'):
    if fn == '':
        continue

    (dir_name,file_name) = path.split(fn)
    (base_name,extension) = path.splitext(file_name)

    new_name = base_name
    write(1,path.join(dir_name,new_name + extension))
