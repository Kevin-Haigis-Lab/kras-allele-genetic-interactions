#!/bin/bash

# Make a new page for my lab notebook

todaysDate=$(date +%Y-%m-%d)
newFileName="./reports/content/post/$todaysDate-$todaysDate.md"

templateFile="---
title: '"$todaysDate"'
author: ~
date: '"$todaysDate"'
slug: '"$todaysDate"'
categories: []
tags: []
subtitle: ''
summary: ''
authors: []
lastmod: '"$todaysDate"'
featured: no
image:
  caption: ''
  focal_point: ''
  preview_only: no
projects: []
---

"

echo "$templateFile" > $newFileName

echo "New notebook page at: $newFileName"
