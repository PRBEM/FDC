#!/bin/bash
for i in {0..23}
do
	sed -i 's/Div7/Div8/g' params_MLT$i.txt
done
