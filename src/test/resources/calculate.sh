#!/bin/bash

cd $1
find . -type f -exec md5sum {} +