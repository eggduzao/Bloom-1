#!/usr/bin/env bash
find . -name "*.DS_Store*" -type f -delete
find . -name "*.Rout" -type f -delete
find . -name "*.RData" -type f -delete
find . -name "*.pyc" -type f -delete
git config --global push.default simple
git config --global user.name "eggduzao"
git config --global user.email eggduzao@gmail.com
git config --global credential.helper 'cache --timeout=3600000000000'
