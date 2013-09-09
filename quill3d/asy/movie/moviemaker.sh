#!/bin/bash

convert -density 100 -delay 15 -loop 0 frame*.png movie-"$1".gif # delay in 1/100 s
