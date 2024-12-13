#!/bin/bash
grep "KRAB" $1 | awk '{if($8 > 13){print}}' 
grep "SCAN" $1 | awk '{if($8 > 10){print}}' 
grep "DUF3669" $1 | awk '{if($8 > 21){print}}' 
grep -v "KRAB" $1 | grep -v "SCAN" | grep -v "DUF3669" | awk '{if($8 > 0){print}}' 
 