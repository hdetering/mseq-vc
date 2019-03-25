#!/bin/awk
# vim: syntax=awk tabstop=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Reformat VarScan VCF output for downstream analysis.
# Formatting steps performed:
#  1. Move INFO:SS to FORMAT:SS
#  2. Move INFO:SSC to FORMAT:SSC
# NOTE: Assuming a two-sample VCF (1st sample healthy. 2nd normal).
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2019-03-08
#------------------------------------------------------------------------------

BEGIN {
  OFS = "\t"; # output fields seperated by TABs
}
# HEADER: INFO field "SS"
# - make it a FORMAT field
# - change Type to "Integer"
/^##INFO=<ID=SS,/ {
  s1 = gensub("INFO", "FORMAT", 1, $0);
  s2 = gensub("Type=String", "Type=Integer", 1, s1);
  $0 = s2;
}
# HEADER: INFO field "SSC"
# - make it a FORMAT field
/^##INFO=<ID=SSC,/ {
  s = gensub("INFO", "FORMAT", 1, $0);
  $0 = s;
}
# print header rows (modified or not)
/^#/ { print; }
# RECORD:
# - move SS, SSC to FORMAT
!/^#/ {
  new_infos = "";
  add_fmt_key = "";
  split("", add_fmt_val); # empty value array
  
  # parse INFO fields
  split($8, infos, ";");
  for (i=1; i<=length(infos); i++) {
    tag = infos[i];
    split(tag, kv, "=");
    if (kv[1] == "SS") { # store value of SS field
      ss = kv[2];
      add_fmt_key = add_fmt_key":"kv[1];
    } else if (kv[1] == "SSC") { # store value of SSC field
      ssc = kv[2];
      add_fmt_key = add_fmt_key":"kv[1];
    } else { # any additional field gets emitted as-is
      if (length(new_infos) == 0) {
        new_infos = tag;
      } else {
        new_infos = new_infos";"tag;
      }
    }
  }
  # replace INFO tags with updated ones
  $8 = new_infos;
  # update FORMAT fields
  split($9, old_fmt_key, ":");
  $9 = $9add_fmt_key
  split($9, new_fmt_key, ":");

  # parse samples' FORMAT values
  # NOTE: assuming 1st sample is normal, 2nd tumor
  split($10, fmt_val_n, ":");
  split($11, fmt_val_t, ":");
  for (i=1; i<=length(old_fmt_key); i++) {
    fmt_n[old_fmt_key[i]] = fmt_val_n[i];
    fmt_t[old_fmt_key[i]] = fmt_val_t[i];
  }
  
  # add FORMAT field values to samples
  # SS (somatic status) field
  if (fmt_n["GT"] == "0/1" || fmt_n["GT"] == "1/1") { # normal is mutated
    fmt_n["SS"] = 1; # germline mutation
  } else { # normal is wildtype (not agreeing with SS=1)
    fmt_n["SS"] = 0; # wildtype
  }
  if (fmt_t["GT"] == "0/1" || fmt_n["GT"] == "1/1") { # tumor is mutated
    fmt_t["SS"] = ss; # original SS should have been 1(germline) or 2(somatic)
  } else {
    if (fmt_n["SS"] == 1) {
      fmt_t["SS"] = 3; # LOH
    } else {
      fmt_t["SS"] = 0; # wildtype
    }
  }
  # SSC (somatic score) field
  fmt_n["SSC"] = ssc;
  fmt_t["SSC"] = ssc;
  # concatenate new set of FORMAT fields
  fmt_str_n = fmt_n[new_fmt_key[1]];
  fmt_str_t = fmt_t[new_fmt_key[1]];
  for (i=2; i<=length(new_fmt_key); i++) {
    fmt_str_n = fmt_str_n":"fmt_n[new_fmt_key[i]];
    fmt_str_t = fmt_str_t":"fmt_t[new_fmt_key[i]];
  }
  $10 = fmt_str_n;
  $11 = fmt_str_t;

  print;
}
