sed -E 's/dimension: 3/dimension: 4/;s/(sizes:.*)/\1 3/;s/(space directions: .*)/\1 none/;s/(kinds:.*)/\1 RGB-color/;s/data file: (.*)red(.*)/data file: LIST 3\n\1red\2\n\1green\2\n\1blue\2/' $x
