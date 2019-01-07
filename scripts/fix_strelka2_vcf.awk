#!/bin/awk
# Replaces empty columns with '.'
BEGIN {
  FS = "\t"
  OFS = "\t"
}
!/^#/ { # fix non-header lines
  for (i=1; i<=NF; i++) {
    if (length($i) == 0) {
      $i = "."
    }
  }
  print
  next
}
{ # print all other lines (headers etc.)
  print 
}
