#!/bin/bash

tmpdir=$(mktemp -d)
wget -nv -P "$tmpdir" -r -l 1 http://www.win.tue.nl/~aeb/graphs/srg/srgtab.html
find $tmpdir -iname "srgtab?*html" | while read filename; do
cat "$filename" | grep -e "^<tr bgcolor" -A 1 |
    while read l; do
	read comment
	echo "$l $comment "|
	    perl -ape "s/.*bgcolor=.(.......)....../\1|/g" |
	    perl -ape "s/<\/td> <td>/|/g" |
	    sed "s/&nbsp.//g" |
	    sed "s/<\/td>//g" |
	    sed "s/<\/tr>//g"
    done
done  > brouwer.tmp
