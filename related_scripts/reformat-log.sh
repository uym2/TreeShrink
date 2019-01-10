#!/bin/bash

cat $1|grep "will be cut" -B1|tr '\n' '@'|sed -e "s/:@//g"|tr '@' '\n'|awk '{print $1,$6,$13}'
