#!/bin/csh -f

# compute complexity of 3D object according to how similar the 2D views are
# used for my SHREC21 submission

# NOTE: the models in the SHREC21 datasets were sometimes numbered from 0 upwards or from 1
# this requires manually adjusting the initial value of c below

# delete any old results
if -f complexityView1.dat then
	rm complexityView1.dat
endif
if -e tmp/complexityView1.dat then
	rm tmp/complexityView1.dat
endif

foreach i (*pix)
	./ramer -i $i -o ./tmp/lines -d 0.5
	./super2pixel -i ./tmp/lines -o {$i:r}.pixRamer
	./pix_transf -i {$i:r}.pixRamer -o {$i:r}.pixRamerBasic -b
end

# count number of models
set count = 0
foreach i (*_030.pix)
	@ count++
end

echo > blank.txt

# if models start at 1
set c = 1
# if models start at 0 (otherwise comment out the following line)
set c = 0
while ($c < $count)
	#echo ===== $c ====
	if -e ./tmp/scores then
		rm ./tmp/scores
	endif
	foreach i ({$c}_*.pixRamerBasic 0{$c}_*.pixRamerBasic 00{$c}_*.pixRamerBasic)
		foreach j ({$c}_*.pixRamerBasic 0{$c}_*.pixRamerBasic 00{$c}_*.pixRamerBasic)
			if ($i != $j) then
				#echo $i $j
				cat $i blank.txt $j > ./tmp/poly
			
				./match_poly_hausdorf < ./tmp/poly >> ./tmp/scores
			endif
		end
	end
	sed -f sedfile ./tmp/scores > ./tmp/scores2
	./avgPlain < ./tmp/scores2 >> ./tmp/complexityView1.dat

	@ c++
end

./delete_CR2 ./tmp/complexityView1.dat complexityView1.dat
