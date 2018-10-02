#!/usr/bin/env bash
if [[ ($# -ne 1) && ($# -ne 2) ]]
then
    echo "usage is $0 <grep pattern to look for in 'jar tvf' output> [<top-of-dir-tree> or, if missing, current dir]"
else
    THING_TO_LOOKFOR="$1"
    DIR=${2:-.}
    if [ ! -d $DIR ]; then
        echo "directory [$DIR] does not exist";
        exit 1;
    fi
    find "$DIR" -iname \*.jar | while read f ; do (jar tf $f | awk '{print "'"$f"'" "  " $0}' | grep -i "$THING_TO_LOOKFOR") ; done
fi