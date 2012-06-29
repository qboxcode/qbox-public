#!/bin/bash
awk '/kpoint/ {print $3,$4,$5}' $1
